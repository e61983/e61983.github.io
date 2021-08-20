package main

import (
	"fmt"
	"math"
)

const (
	a    float64 = 6378137.0
	b    float64 = 6356752.3142451
	lon0 float64 = 121 * math.Pi / 180
	k0   float64 = 0.9999
	dx   float64 = 250000
	dy   float64 = 0
)

var (
	e  float64 = 1 - math.Pow(b, 2)/math.Pow(a, 2)
	e2 float64 = (1 - math.Pow(b, 2)/math.Pow(a, 2)) / (math.Pow(b, 2) / math.Pow(a, 2))
)

func LonLat2TM2(lon, lat float64) (x, y float64) {

	lon = (lon - math.Floor((lon+180)/360)*360) * math.Pi / 180
	lat = lat * math.Pi / 180

	V := a / math.Sqrt(1-e*math.Pow(math.Sin(lat), 2))
	T := math.Pow(math.Tan(lat), 2)
	C := e2 * math.Pow(math.Cos(lat), 2)
	A := math.Cos(lat) * (lon - lon0)
	M := a * ((1.0-e/4.0-3.0*math.Pow(e, 2)/64.0-5.0*math.Pow(e, 3)/256.0)*lat -
		(3.0*e/8.0+3.0*math.Pow(e, 2)/32.0+45.0*math.Pow(e, 3)/1024.0)*
			math.Sin(2.0*lat) + (15.0*math.Pow(e, 2)/256.0+45.0*math.Pow(e, 3)/1024.0)*
		math.Sin(4.0*lat) - (35.0*math.Pow(e, 3)/3072.0)*math.Sin(6.0*lat))

	x = dx + k0*V*(A+(1-T+C)*math.Pow(A, 3)/6+(5-18*T+math.Pow(T, 2)+72*C-58*e2)*math.Pow(A, 5)/120)

	y = dy + k0*(M+V*math.Tan(lat)*(math.Pow(A, 2)/2+(5-T+9*C+4*math.Pow(C, 2))*math.Pow(A, 4)/24+(61-58*T+math.Pow(T, 2)+600*C-330*e2)*math.Pow(A, 6)/720))
	return
}

func TM22LonLat(x, y float64) (lon, lat float64) {
	x -= dx
	y -= dy

	// Calculate the Meridional Arc
	M := y / k0

	// Calculate Footprint Latitude
	mu := M / (a * (1.0 - e/4.0 - 3*math.Pow(e, 2)/64.0 - 5*math.Pow(e, 3)/256.0))

	e1 := (1.0 - math.Sqrt(1.0-e)) / (1.0 + math.Sqrt(1.0-e))

	J1 := (3*e1/2 - 27*math.Pow(e1, 3)/32.0)
	J2 := (21*math.Pow(e1, 2)/16 - 55*math.Pow(e1, 4)/32.0)
	J3 := (151 * math.Pow(e1, 3) / 96.0)
	J4 := (1097 * math.Pow(e1, 4) / 512.0)

	fp := mu + J1*math.Sin(2*mu) + J2*math.Sin(4*mu) + J3*math.Sin(6*mu) + J4*math.Sin(8*mu)

	// Calculate Latitude and Longitude

	C1 := e2 * math.Pow(math.Cos(fp), 2)
	T1 := math.Pow(math.Tan(fp), 2)
	R1 := a * (1 - e) / math.Pow((1-e*math.Pow(math.Sin(fp), 2)), (3.0/2.0))
	N1 := a / math.Pow((1-e*math.Pow(math.Sin(fp), 2)), 0.5)

	D := x / (N1 * k0)

	// 計算緯度
	Q1 := N1 * math.Tan(fp) / R1
	Q2 := (math.Pow(D, 2) / 2.0)
	Q3 := (5 + 3*T1 + 10*C1 - 4*math.Pow(C1, 2) - 9*e2) * math.Pow(D, 4) / 24.0
	Q4 := (61 + 90*T1 + 298*C1 + 45*math.Pow(T1, 2) - 3*math.Pow(C1, 2) - 252*e2) * math.Pow(D, 6) / 720.0
	lat = fp - Q1*(Q2-Q3+Q4)

	// 計算經度
	Q5 := D
	Q6 := (1 + 2*T1 + C1) * math.Pow(D, 3) / 6
	Q7 := (5 - 2*C1 + 28*T1 - 3*math.Pow(C1, 2) + 8*e2 + 24*math.Pow(T1, 2)) * math.Pow(D, 5) / 120.0
	lon = lon0 + (Q5-Q6+Q7)/math.Cos(fp)

	lat = (lat * 180) / math.Pi //緯
	lon = (lon * 180) / math.Pi //經

	return
}

func TWD672TWD97(x, y float64) (x_97, y_97 float64) {
	const A float64 = 0.00001549
	const B float64 = 0.000006521
	x_97 = x + 807.8 + A*x + B*y
	y_97 = y - 248.6 + A*y + B*x
	return
}

func TWD972TWD67(x, y float64) (x_67, y_67 float64) {
	const A float64 = 0.00001549
	const B float64 = 0.000006521
	x_67 = x - 807.8 - A*x - B*y
	y_67 = y + 248.6 - A*y - B*x
	return
}

// References:
// https://www.sunriver.com.tw/taiwanmap/grid_tm2_convert.php
func main() {
	const x_67 float64 = 247342
	const y_67 float64 = 2652336

	x_97, y_97 := TWD672TWD97(x_67, y_67)

	fmt.Printf("TWD67:\n")
	fmt.Printf("\t%f, %f\n", x_67, y_67)
	fmt.Printf("TWD97:\n")
	fmt.Printf("\t%f, %f\n", x_97, y_97)

	lon, lat := TM22LonLat(x_97, y_97)

	fmt.Printf("LonLat:\n")
	fmt.Printf("\t%f, %f\n", lon, lat)
}
