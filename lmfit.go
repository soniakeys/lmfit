// Public domain.

// Package lmfit fits sky observations to linear motion along a great circle.
package lmfit

import (
	"math"

	"github.com/soniakeys/coord"
	"github.com/soniakeys/unit"
)

/*
const (
	twoPi     = 2 * math.Pi
	threePi   = 3 * math.Pi
	arcSecRad = math.Pi / (180 * 3600)
)
*/
// LmFit represents the linear motion resulting from a great circle fit
// of a sequence of spherical coordinates.
// It can be queried for various statistics on the fit and used to generate
// synthetic coordinates corresponding to the linear motion.
type LmFit struct {
	// various intermediate values, useful for generating statistics after
	// the fit is done
	mRot  coord.M3    // rotation matrix
	lon0  unit.Angle  // lon (ra) offset
	rs    coord.SphrS // rotated ra and dec
	t0    float64     // time offset
	nTime []float64   // normalized times
	// fit solution parameters
	r0, d0 unit.Angle
	rr, dr unit.Angle // actually angle/day
}

// New does the great circle fit.
//
// Args:
//   mjd  -- mjd times of observations
//   s    -- ra, dec coordinates
func New(mjd []float64, e coord.EquaS) *LmFit {
	nObs := len(mjd)
	if nObs < 2 || len(e) != nObs {
		return nil
	}
	var lmf LmFit

	// convert obs to cartesian
	var c coord.CartS
	c.FromEquaS(e)

	// vector normal to motion
	var n coord.Cart
	n.Cross(&c[0], &c[nObs-1])
	nmag2 := n.Square()
	nmag := math.Sqrt(nmag2)

	// rotation angle is angle from n to +z
	xy := math.Hypot(n.X, n.Y)

	if tana := xy / n.Z; !(math.Abs(tana) > .0003) {
		// if n is close to the pole already, don't mess with rotation.
		lmf.mRot = coord.M3{1, 0, 0, 0, 1, 0, 0, 0, 1} // identity
		lmf.rs = e.SphrS()                             // copy
	} else {
		// rotation axis is right angle to n projected onto z=0 plane
		// and normalized to unit vector
		sina := xy / nmag
		cosa := n.Z / nmag
		gcix := n.Y / xy
		gciy := -n.X / xy

		// build mRot, rotation matrix that will rotate the coordinate system
		// to the obs, so that a cylindrical projection will have negligible
		// distortion.)
		sinagx := sina * gcix
		sinagy := sina * gciy
		onemcosa := 1 - cosa
		onemcosagx := onemcosa * gcix
		onemcosagxgy := onemcosagx * gciy
		lmf.mRot[0] = cosa + onemcosagx*gcix
		lmf.mRot[1] = onemcosagxgy
		lmf.mRot[2] = sinagy
		lmf.mRot[3] = onemcosagxgy
		lmf.mRot[4] = cosa + onemcosa*gciy*gciy
		lmf.mRot[5] = -sinagx
		lmf.mRot[6] = -sinagy
		lmf.mRot[7] = sinagx
		lmf.mRot[8] = cosa

		// rotate all of cart
		var rotated coord.CartS
		rotated.Mult3S(&lmf.mRot, c)

		// transpose rotation array so it will derotate, after
		// least-squares fit.
		lmf.mRot.Transpose(&lmf.mRot)
		/*
			// invariant:  first and last points should be on the z=0 plane
			// of the rotated coordinate system now.
			const rTol = arcSecRad
			if math.Abs(rotated[0].Z) > rTol || math.Abs(rotated[nObs-1].Z) > rTol {
				panic(fmt.Sprintf("%f %f %f\n*** rotation failed ***",
					rotated[0].Z, rotated[nObs-1].Z, rTol))
			}
		*/
		// convert back to spherical coordinates for adjustment.
		// (this does the cylindrical projection.)
		lmf.rs.FromCartS(rotated)
	}

	// normalize longitude (ra) to near 0 to avoid wraparound problems
	lon0 := lmf.rs[0].Lon
	lmf.lon0 = lon0
	for i, r1 := range lmf.rs {
		lmf.rs[i].Lon = (r1.Lon - lon0 + math.Pi).Mod1() - math.Pi
	}

	// normalize time to near 0 to maintain precision
	t0 := mjd[0]
	lmf.t0 = t0
	lmf.nTime = make([]float64, nObs)
	for i, m1 := range mjd {
		lmf.nTime[i] = m1 - t0
	}

	if nObs == 2 {
		lmf.r0 = 0
		lmf.rr = lmf.rs[1].Lon.Div(lmf.nTime[1])
		lmf.d0 = 0
		lmf.dr = 0
	} else {
		// here's the least squares stuff
		var sumT, sumT2 float64
		var sumLon, sumLat, sumTLon, sumTLat unit.Angle
		for i, t1 := range lmf.nTime {
			sumT += t1
			sumLon += lmf.rs[i].Lon
			sumLat += lmf.rs[i].Lat
			sumT2 += t1 * t1
			sumTLon += lmf.rs[i].Lon.Mul(t1)
			sumTLat += lmf.rs[i].Lat.Mul(t1)
		}
		fn := float64(nObs)
		d := fn*sumT2 - sumT*sumT
		lmf.r0 = (sumLon.Mul(sumT2) - sumTLon.Mul(sumT)).Div(d)
		lmf.rr = (sumTLon.Mul(fn) - sumLon.Mul(sumT)).Div(d)
		lmf.d0 = (sumLat.Mul(sumT2) - sumTLat.Mul(sumT)).Div(d)
		lmf.dr = (sumTLat.Mul(fn) - sumLat.Mul(sumT)).Div(d)
	}
	return &lmf
}

// Pos returns position at time t with fitted linear motion.
func (lmf *LmFit) Pos(t float64) *coord.Equa {
	nt := t - lmf.t0
	s := &coord.Sphr{
		Lon: lmf.r0 + lmf.rr.Mul(nt) + lmf.lon0,
		Lat: lmf.d0 + lmf.dr.Mul(nt),
	}
	var c coord.Cart
	var p coord.Equa
	return p.FromCart(c.Mult3(&lmf.mRot, c.FromSphr(s)))
}

// Res computes and returns residuals, as documented at RmsRes.
func (lmf *LmFit) Res() coord.SphrS {
	// residuals are computed on rotated-and-derotated observed to
	// minimize any systematic errors from rotating and derotating the
	// computed.

	// observed positions:  copy from rotated observed and fix up ra
	so := append(coord.SphrS{}, lmf.rs...)
	for i := range so {
		so[i].Lon += lmf.lon0
	}
	// then derotate back up to original place in the sky
	var c coord.CartS
	so.FromCartS(c.Mult3S(&lmf.mRot, c.FromSphrS(so)))

	// computed positions, derotated similarly
	sc := make(coord.SphrS, len(so))
	for i, nt := range lmf.nTime {
		sc[i] = coord.Sphr{
			Lon: lmf.r0 + lmf.rr.Mul(nt) + lmf.lon0,
			Lat: lmf.d0 + lmf.dr.Mul(nt),
		}
	}
	sc.FromCartS(c.Mult3S(&lmf.mRot, c.FromSphrS(sc)))

	// repurpose so to hold residuals
	for i, so1 := range so {
		so[i] = coord.Sphr{
			Lon: (so1.Lon - sc[i].Lon).Mul(sc[i].Lat.Cos()),
			Lat: so1.Lat - sc[i].Lat,
		}
	}
	return so
}

// RmsRes returns rms of 2d residuals and the component Lon and dec residuals.
//
// Note:  The 2d residual is the 2d distance between computed and observed.
// By contrast, it is not the rms of the ra and dec residuals considered as
// separate values.  (The 2d distance that this function returns is smaller
// than the rms of separate values by a factor of sqrt(2).)
func (lmf *LmFit) RmsRes() (unit.Angle, coord.SphrS) {
	res := lmf.Res()
	var s float64
	for _, r1 := range res {
		s += (r1.Lon*r1.Lon + r1.Lat*r1.Lat).Rad()
	}
	return unit.Angle(math.Sqrt(s / float64(len(res)))), res
}

// Rms returns just the rms, as documented at RmsRes.
func (lmf *LmFit) Rms() unit.Angle {
	rms, _ := lmf.RmsRes()
	return rms
}
