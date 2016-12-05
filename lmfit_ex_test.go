// Public domain.

// Example shows New and LmFit.Pos
package lmfit_test

import (
	"fmt"

	"github.com/soniakeys/coord"
	"github.com/soniakeys/lmfit"
	"github.com/soniakeys/sexagesimal"
	"github.com/soniakeys/unit"
)

func equaString(prec int, s *coord.Equa) string {
	return fmt.Sprintf("{RA %2.*s, Dec %2.*s}",
		prec, sexa.FmtRA(s.RA),
		prec, sexa.FmtAngle(s.Dec))
}

func Example() {
	// New, with two observations
	mjd := []float64{56123, 56123.01}
	s := coord.EquaS{
		{RA: 0, Dec: unit.NewAngle(' ', 89, 59, 40)},
		{RA: 0, Dec: unit.NewAngle(' ', 90, 0, 0)},
	}
	f := lmfit.New(mjd, s)

	// Reproduce positions at first, last times.
	// Between, a synthesized position for the mean time.
	fmt.Println("degrees:")
	fmt.Println(equaString(1, f.Pos(mjd[0])))
	fmt.Println(equaString(1, f.Pos((mjd[0]+mjd[1])/2)))
	fmt.Println(equaString(1, f.Pos(mjd[1])))
	// Output:
	// degrees:
	// {RA   0ʰ 0ᵐ 0.0ˢ, Dec  89°59′40.0″}
	// {RA   0ʰ 0ᵐ 0.0ˢ, Dec  89°59′50.0″}
	// {RA   0ʰ 0ᵐ 0.0ˢ, Dec  90° 0′ 0.0″}
}
