// Public domain.

package lmfit_test

import (
	"fmt"
	"math"

	"github.com/soniakeys/coord"
	"github.com/soniakeys/lmfit"
	"github.com/soniakeys/unit"
)

func ExampleLmFit_RmsRes() {
	// Three observations, the middle one an arc second off in both RA and Dec.
	mjd := []float64{56123, 56123.01, 56123.02}
	s := coord.EquaS{
		{
			RA:  unit.NewRA(23, 59, 59),
			Dec: 0,
		},
		{
			RA:  unit.NewRA(0, 0, 1./15),
			Dec: unit.NewAngle(' ', 0, 0, 1),
		},
		{
			RA:  unit.NewRA(0, 0, 1),
			Dec: 0,
		},
	}
	f := lmfit.New(mjd, s)

	// show fitted positions
	fmt.Println("fitted positions, in degrees:")
	for _, t := range mjd {
		fmt.Println(equaString(2, f.Pos(t)))
	}

	// Get both RMS and residuals
	rms, res := f.RmsRes()

	fmt.Println("\nresiduals, in arc seconds:")
	for _, r := range res {
		fmt.Printf("{RA %.2f, Dec %.2f}\n", r.Lon.Sec(), r.Lat.Sec())
	}
	fmt.Printf("\nrms: %.2f\n", rms.Sec())

	// compute RMS by hand to illustrate formula
	ot := 1. / 3                 // one third (of an arc second)
	ots := ot * ot               // one third squared
	res1 := math.Sqrt(ots + ots) // 2d residual of first point
	tt := 2. / 3
	tts := tt * tt
	res2 := math.Sqrt(tts + tts) // of second point
	res3 := res1                 // of third
	fmt.Printf("rms: %.2f\n", math.Sqrt((res1*res1+res2*res2+res3*res3)/3))
	// this is equivalent:
	fmt.Printf("rms: %.2f\n", math.Sqrt((ots+ots+tts+tts+ots+ots)/3))
	// this is not:
	fmt.Printf("not: %.2f\n", math.Sqrt((ots+ots+tts+tts+ots+ots)/6))
	// Output:
	// fitted positions, in degrees:
	// {RA  23ʰ59ᵐ59.02ˢ, Dec   0° 0′ 0.33″}
	// {RA   0ʰ 0ᵐ 0.02ˢ, Dec   0° 0′ 0.33″}
	// {RA   0ʰ 0ᵐ 1.02ˢ, Dec   0° 0′ 0.33″}
	//
	// residuals, in arc seconds:
	// {RA -0.33, Dec -0.33}
	// {RA 0.67, Dec 0.67}
	// {RA -0.33, Dec -0.33}
	//
	// rms: 0.67
	// rms: 0.67
	// rms: 0.67
	// not: 0.47
}
