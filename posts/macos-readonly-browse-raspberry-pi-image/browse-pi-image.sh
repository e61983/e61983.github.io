#!/usr/bin/env bash
#
set -euo pipefail

# ── Constants ────────────────────────────────────────
readonly DEFAULT_IMAGE_NAME="image"
readonly CONTAINER_IMAGE="debian:stable-slim"

# ── Utility Functions ────────────────────────────────
script_dir() {
  cd "$(dirname "${BASH_SOURCE[0]}")" && pwd
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [--selftest] [image-path]

  Browse a Raspberry Pi disk image root filesystem in read-only mode.
  If no path is provided, defaults to ${DEFAULT_IMAGE_NAME} in the script directory.

Options:
  --selftest   Non-interactive self-test: mount and verify rootfs content and read-only mode, then exit.
  -h, --help   Show this help message.

Examples:
  $(basename "$0")
  $(basename "$0") ~/Downloads/other.img
EOF
}

# Verify image file exists and Docker daemon is running
preflight() {
  local image="$1"
  if [[ ! -f "$image" ]]; then
    echo "Error: Image file not found: $image" >&2
    echo >&2
    usage >&2
    return 1
  fi
  if ! docker info >/dev/null 2>&1; then
    echo "Error: Docker daemon is not running. Please start Docker Desktop first." >&2
    return 1
  fi
}

# Generate setup script to run inside container: detect partitions → mount read-only → self-test / interactive shell
container_setup_script() {
  cat <<'SETUP'
set -euo pipefail

LOOP=""
cleanup() {
  [ -n "$LOOP" ] || return 0
  cd / 2>/dev/null || true
  umount -R /rootfs 2>/dev/null || true
  losetup -d "$LOOP" 2>/dev/null || true
}
trap cleanup EXIT

echo "[setup] Creating read-only loop device and detecting partitions..."
LOOP=$(losetup -fP --read-only --show /image.img)
partx -a "$LOOP" 2>/dev/null || true

# When udev is not running in container, partx cannot create device nodes → create manually from /proc/partitions
for _p in 1 2; do
  _dev="${LOOP}p${_p}"
  if [ ! -e "$_dev" ]; then
    _name="${LOOP#/dev/}p${_p}"
    _maj=$(awk -v n="$_name" '$4==n{print $1}' /proc/partitions)
    _min=$(awk -v n="$_name" '$4==n{print $2}' /proc/partitions)
    [ -n "$_maj" ] && mknod "$_dev" b "$_maj" "$_min" 2>/dev/null || true
  fi
done

ROOT_DEV="${LOOP}p2"
BOOT_DEV="${LOOP}p1"

if [ ! -e "$ROOT_DEV" ]; then
  echo "[setup] Rootfs partition not found ($ROOT_DEV). Diagnostic info:" >&2
  lsblk "$LOOP" 2>/dev/null >&2 || losetup -l >&2 || true
  exit 1
fi

mkdir -p /rootfs
echo "[setup] Mounting rootfs read-only ($ROOT_DEV)..."
if ! mount -o ro "$ROOT_DEV" /rootfs 2>/dev/null; then
  echo "[setup] Direct mount failed, trying with ro,noload to skip journal replay..." >&2
  if ! mount -o ro,noload "$ROOT_DEV" /rootfs; then
    echo "[setup] Failed to mount rootfs. Diagnostic info:" >&2
    lsblk "$LOOP" 2>/dev/null >&2 || true
    exit 1
  fi
fi

# Mount boot partition at rootfs boot mount point
if [ -d /rootfs/boot/firmware ]; then
  mount -o ro "$BOOT_DEV" /rootfs/boot/firmware 2>/dev/null || true
elif [ -d /rootfs/boot ]; then
  mount -o ro "$BOOT_DEV" /rootfs/boot 2>/dev/null || true
fi

# ── Self-test mode ──
if [ -n "${SELFTEST:-}" ]; then
  echo "[selftest] --- ls /rootfs ---"
  ls -1 /rootfs
  echo "[selftest] --- /etc/os-release ---"
  cat /rootfs/etc/os-release 2>/dev/null || { echo "FAIL: /etc/os-release not found"; exit 1; }
  echo "[selftest] --- Read-only check ---"
  if touch /rootfs/__rw_test__ 2>/dev/null; then
    echo "FAIL: rootfs is writable!"; rm -f /rootfs/__rw_test__ 2>/dev/null || true; exit 1
  fi
  echo "OK: rootfs is read-only"
  echo "[selftest] --- Export check ---"
  cp /rootfs/etc/hostname /export/__selftest_hostname__
  echo "OK: Exported to /export"
  echo "[selftest] All checks passed"
  exit 0
fi

# ── Interactive browse mode ──
cat <<'BANNER'

────────────────────────────────────────────────────────────
 Raspberry Pi image — Read-only browse shell
────────────────────────────────────────────────────────────
 • rootfs mounted at /rootfs   (read-only throughout, original image not modified)
 • Copy back to Mac:   cp /rootfs/<path> /export/
 • Exit:       exit   or   Ctrl-D
────────────────────────────────────────────────────────────

BANNER

cat > /root/.browserc <<'RC'
cd /rootfs
PS1='pi-image:\w\$ '
RC
bash --rcfile /root/.browserc -i || true
SETUP
}

# Launch privileged container to mount image read-only and run setup script
run_container() {
  local image="$1" export_dir="$2" selftest="$3"
  mkdir -p "$export_dir"

  local setup_script
  setup_script="$(container_setup_script)"

  local -a flags=( --rm --privileged
    -v "${image}:/image.img:ro"
    -v "${export_dir}:/export:rw" )

  if [[ "$selftest" == "1" ]]; then
    flags+=( -e SELFTEST=1 )
  else
    flags+=( -it )
  fi

  docker run "${flags[@]}" "$CONTAINER_IMAGE" bash -c "$setup_script"
}

# ── Main routine ─────────────────────────────────────
main() {
  local selftest=0
  local -a positional=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -h|--help) usage; exit 0 ;;
      --selftest) selftest=1; shift ;;
      --) shift; while [[ $# -gt 0 ]]; do positional+=("$1"); shift; done; break ;;
      -*) echo "Unknown option: $1" >&2; echo >&2; usage >&2; exit 1 ;;
      *) positional+=("$1"); shift ;;
    esac
  done

  local image
  if [[ ${#positional[@]} -gt 0 ]]; then
    image="${positional[0]}"
  else
    image="$(script_dir)/${DEFAULT_IMAGE_NAME}"
  fi

  preflight "$image"

  # docker -v requires absolute path
  image="$(cd "$(dirname "$image")" && pwd)/$(basename "$image")"

  local export_dir
  export_dir="$(script_dir)/_export"

  run_container "$image" "$export_dir" "$selftest"
}

main "$@"
