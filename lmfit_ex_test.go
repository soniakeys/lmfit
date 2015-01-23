// Public domain.

// Example shows New and LmFit.Pos
package lmfit_test

import (
	"fmt"

	"github.com/soniakeys/coord"
	"github.com/soniakeys/lmfit"
	"github.com/soniakeys/sexagesimal"
)

func sphrString(prec int, s *coord.Sphr) string {
	return fmt.Sprintf("{RA %2.*s, Dec %2.*s}",
		prec, sexa.NewFmtRA(s.RA),
		prec, sexa.NewFmtAngle(s.Dec))
}

func Example() {
	// New, with two observations
	mjd := []float64{56123, 56123.01}
	s := coord.SphrS{
		{RA: 0, Dec: sexa.NewAngle(false, 89, 59, 40).Rad()},
		{RA: 0, Dec: sexa.NewAngle(false, 90, 0, 0).Rad()},
	}
	f := lmfit.New(mjd, s)

	// Reproduce positions at first, last times.
	// Between, a synthesized position for the mean time.
	fmt.Println("degrees:")
	fmt.Println(sphrString(1, f.Pos(mjd[0])))
	fmt.Println(sphrString(1, f.Pos((mjd[0]+mjd[1])/2)))
	fmt.Println(sphrString(1, f.Pos(mjd[1])))
	// Output:
	// degrees:
	// {RA  0ʰ 0ᵐ 0.0ˢ, Dec 89°59′40.0″}
	// {RA  0ʰ 0ᵐ 0.0ˢ, Dec 89°59′50.0″}
	// {RA  0ʰ 0ᵐ 0.0ˢ, Dec 90° 0′ 0.0″}
}
