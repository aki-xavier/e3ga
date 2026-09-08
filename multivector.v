module e3ga

// Plain 3D Euclidean Geometric Algebra Cl(3,0) — the "ordinary" GA, in
// contrast to the 5D conformal algebra of the parent cga module.
//
// A multivector has 8 components.  Each basis blade is indexed by its mask
// (bit 0 = e1, bit 1 = e2, bit 2 = e3), so the geometric product of two basis
// blades is a pure XOR/mask-sign computation:
//
//   [0] grade 0: 1
//   [1] [2] [4]  grade 1: e1, e2, e3
//   [3] [5] [6]  grade 2: e12, e13, e23
//   [7] grade 3: I = e123   (pseudoscalar, I^2 = -1)
//
// Metric: e1^2 = e2^2 = e3^2 = +1 and ei ej = -ej ei for i != j.  The even
// subalgebra (scalar + bivector) is isomorphic to the quaternions: a unit
// rotor R = cos(t/2) - sin(t/2) . Bhat rotates by v' = R v R~ (apply).
import math

pub const num_components = 8
pub const num_grades = 4

pub struct Multivector {
pub mut:
	values [8]f64
}

// --- constructors -----------------------------------------------------------

// mv_zero returns the zero multivector.
pub fn mv_zero() Multivector {
	return Multivector{}
}

// mv_scalar returns a scalar multivector (only the grade-0 component).
pub fn mv_scalar(s f64) Multivector {
	mut m := Multivector{}
	m.values[0] = s
	return m
}

// mv_vector builds a vector (grade 1) from euclidean components (x, y, z).
pub fn mv_vector(x f64, y f64, z f64) Multivector {
	mut m := Multivector{}
	m.values[1] = x
	m.values[2] = y
	m.values[4] = z
	return m
}

// mv_bivector builds a bivector (grade 2) from canonical components
// (e12, e13, e23).
pub fn mv_bivector(b12 f64, b13 f64, b23 f64) Multivector {
	mut m := Multivector{}
	m.values[3] = b12
	m.values[5] = b13
	m.values[6] = b23
	return m
}

// pseudoscalar returns I = e1 e2 e3 (squares to -1).
pub fn pseudoscalar() Multivector {
	mut m := Multivector{}
	m.values[7] = 1.0
	return m
}

// Basis vectors (grade 1).  Fresh values per call.
pub fn e1() Multivector {
	return mv_vector(1.0, 0.0, 0.0)
}

pub fn e2() Multivector {
	return mv_vector(0.0, 1.0, 0.0)
}

pub fn e3() Multivector {
	return mv_vector(0.0, 0.0, 1.0)
}

// --- component access -------------------------------------------------------

// grade returns the grade-g projection.
pub fn (m Multivector) grade(g int) Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		if popcount(i) == g {
			res.values[i] = m.values[i]
		}
	}
	return res
}

// scalar_part returns the grade-0 component.
pub fn (m Multivector) scalar_part() f64 {
	return m.values[0]
}

// vector_part returns the three grade-1 components (x, y, z).
pub fn (m Multivector) vector_part() [3]f64 {
	return [m.values[1], m.values[2], m.values[4]]!
}

// bivector_part returns the three grade-2 components (e12, e13, e23).
pub fn (m Multivector) bivector_part() [3]f64 {
	return [m.values[3], m.values[5], m.values[6]]!
}

// pseudoscalar_part returns the grade-3 component.
pub fn (m Multivector) pseudoscalar_part() f64 {
	return m.values[7]
}

// is_zero reports whether all components are approximately zero.
pub fn (m Multivector) is_zero() bool {
	for v in m.values {
		if math.abs(v) > 1e-10 {
			return false
		}
	}
	return true
}

// vmax returns the maximum absolute component.
pub fn (m Multivector) vmax() f64 {
	mut mx := 0.0
	for v in m.values {
		a := math.abs(v)
		if a > mx {
			mx = a
		}
	}
	return mx
}

// --- operators --------------------------------------------------------------

pub fn (m Multivector) add(o Multivector) Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		res.values[i] = m.values[i] + o.values[i]
	}
	return res
}

pub fn (m Multivector) sub(o Multivector) Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		res.values[i] = m.values[i] - o.values[i]
	}
	return res
}

pub fn (m Multivector) mul_scalar(s f64) Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		res.values[i] = m.values[i] * s
	}
	return res
}

pub fn (m Multivector) div_scalar(s f64) Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		res.values[i] = m.values[i] / s
	}
	return res
}

pub fn (m Multivector) neg() Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		res.values[i] = -m.values[i]
	}
	return res
}

// eq performs approximate equality (allclose, atol=1e-6).
pub fn (m Multivector) eq(o Multivector) bool {
	for i in 0 .. num_components {
		if math.abs(m.values[i] - o.values[i]) > 1e-6 {
			return false
		}
	}
	return true
}

// copy returns a fresh multivector with the same components.
pub fn (m Multivector) copy() Multivector {
	return m
}

// str renders non-zero components by grade.
pub fn (m Multivector) str() string {
	mut parts := []string{}
	for g in 0 .. num_grades {
		for i in 0 .. num_components {
			if popcount(i) != g {
				continue
			}
			v := m.values[i]
			if math.abs(v) <= 1e-10 {
				continue
			}
			name := blade_name(i)
			if name == '1' {
				parts << '${v:.4f}'
			} else {
				parts << '${v:+.4f}*${name}'
			}
		}
	}
	if parts.len == 0 {
		return 'Multivector(0)'
	}
	return 'Multivector(${parts.join(' ')})'
}

// --- algebra ----------------------------------------------------------------

// gp computes the geometric product.  Basis-blade masks multiply by XOR; the
// sign counts inversions: each pair (i in b, j in a) with j > i contributes a
// swap.  Common basis factors square to +1 (positive-definite metric) and are
// removed by the XOR.
pub fn (m Multivector) gp(o Multivector) Multivector {
	mut res := Multivector{}
	for ma in 0 .. num_components {
		a := m.values[ma]
		if a == 0.0 {
			continue
		}
		for mb in 0 .. num_components {
			b := o.values[mb]
			if b == 0.0 {
				continue
			}
			dst, sign := gp_blade(ma, mb)
			res.values[dst] += sign * a * b
		}
	}
	return res
}

// gp_blade returns the product of two basis blades as (result mask, scalar sign).
fn gp_blade(ma int, mb int) (int, f64) {
	mut swaps := 0
	for i in 0 .. 3 {
		if mb & (1 << i) != 0 {
			swaps += popcount(ma >> (i + 1))
		}
	}
	mut sign := 1.0
	if swaps % 2 != 0 {
		sign = -1.0
	}
	return ma ^ mb, sign
}

// op computes the outer product a ^ b: zero for blades sharing a basis
// factor, else +/- the union blade (same sign rule as gp).
pub fn (m Multivector) op(o Multivector) Multivector {
	mut res := Multivector{}
	for ma in 0 .. num_components {
		a := m.values[ma]
		if a == 0.0 {
			continue
		}
		for mb in 0 .. num_components {
			b := o.values[mb]
			if b == 0.0 {
				continue
			}
			if ma & mb != 0 {
				continue
			}
			dst, sign := gp_blade(ma, mb)
			res.values[dst] += sign * a * b
		}
	}
	return res
}

// ip computes the Hestenes fat-dot inner product.
//
// A|B = sum_{r,s>=1} < <A>_r <B>_s >_|r-s| ; scalar (grade-0) terms are zero.
pub fn (m Multivector) ip(o Multivector) Multivector {
	mut res := Multivector{}
	for ga in 1 .. num_grades {
		a_g := m.grade(ga)
		if a_g.is_zero() {
			continue
		}
		for gb in 1 .. num_grades {
			b_g := o.grade(gb)
			if b_g.is_zero() {
				continue
			}
			prod := a_g.gp(b_g)
			res = res.add(prod.grade(int(math.abs(f64(gb - ga)))))
		}
	}
	return res
}

// reverse applies the reversal involution: grade-k blade * (-1)^(k(k-1)/2).
pub fn (m Multivector) reverse() Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		k := popcount(i)
		if (k * (k - 1) / 2) % 2 != 0 {
			res.values[i] = -m.values[i]
		} else {
			res.values[i] = m.values[i]
		}
	}
	return res
}

// grade_involution flips the sign of odd-grade components.
pub fn (m Multivector) grade_involution() Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		if popcount(i) % 2 != 0 {
			res.values[i] = -m.values[i]
		} else {
			res.values[i] = m.values[i]
		}
	}
	return res
}

// conjugate is the Clifford conjugate (reverse + grade involution).
pub fn (m Multivector) conjugate() Multivector {
	return m.reverse().grade_involution()
}

// dual is the Hodge dual: multiply by I^-1 = -I (I^2 = -1).
// A plane maps to its normal, e.g. dual(e12) = e3.
pub fn (m Multivector) dual() Multivector {
	return m.gp(pseudoscalar().neg())
}

// undual inverts dual (dual(dual(x)) = -x, so undual = -dual).
pub fn (m Multivector) undual() Multivector {
	return m.dual().neg()
}

// meet returns self v other = (self* ^ other*)*, e.g. meet(e12, e13) = e1.
pub fn (m Multivector) meet(o Multivector) Multivector {
	return m.dual().op(o.dual()).dual()
}

// norm returns the euclidean norm sqrt(|<self.reverse(self)>_0|).
pub fn (m Multivector) norm() f64 {
	s := m.gp(m.reverse()).values[0]
	return math.sqrt(math.abs(s))
}

// normalized returns the unit-norm multivector (or zero if degenerate).
pub fn (m Multivector) normalized() Multivector {
	n := m.norm()
	if n < 1e-12 {
		return mv_zero()
	}
	return m.div_scalar(n)
}

// --- rotors (even subalgebra = quaternions) ---------------------------------

// rotor returns the rotation rotor by `angle` radians about `axis`:
// R = exp(-(angle/2) . B) with B the unit bivector of the axis plane.
// Rotating a vector: v' = R v R~ (apply).
pub fn rotor(axis [3]f64, angle f64) Multivector {
	axis_len := math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2])
	if axis_len < 1e-12 {
		panic('rotor: zero axis')
	}
	half := angle / 2.0
	c := math.cos(half)
	s := math.sin(half)
	// Axis n maps to the plane bivector n I with parts (n3, -n2, n1).
	b12 := axis[2] / axis_len
	b13 := -axis[1] / axis_len
	b23 := axis[0] / axis_len
	return mv_scalar(c).sub(mv_bivector(s * b12, s * b13, s * b23))
}

// apply returns r v r~ for a unit rotor; rotates a vector (or blade) by r.
pub fn (r Multivector) apply(v Multivector) Multivector {
	return r.gp(v).gp(r.reverse())
}

// exp exponentiates an even multivector (scalar + bivector): exp(s + B) =
// e^s (cos|B| + B sin|B|/|B|).  Odd-grade components are rejected; any
// bivector B squares to -|B|^2, which is what the closed form uses.
pub fn (m Multivector) exp() Multivector {
	if m.values[1] != 0.0 || m.values[2] != 0.0 || m.values[4] != 0.0 || m.values[7] != 0.0 {
		panic('exp: only even (scalar + bivector) multivectors are supported')
	}
	b := m.bivector_part()
	b_norm := math.sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
	e := math.exp(m.values[0])
	if b_norm < 1e-12 {
		return mv_scalar(e)
	}
	c := e * math.cos(b_norm)
	k := e * math.sin(b_norm) / b_norm
	return mv_scalar(c).add(mv_bivector(b[0] * k, b[1] * k, b[2] * k))
}

// log returns the bivector B with exp(B) = self for a unit rotor.  Its
// magnitude is the half angle: R = cos(t/2) - sin(t/2) Bhat gives
// log(R) = -(t/2) Bhat.
pub fn (m Multivector) log() Multivector {
	s := m.values[0]
	b := m.bivector_part()
	b_norm := math.sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
	if b_norm < 1e-12 {
		if s > 0.0 {
			return mv_zero()
		}
		panic('log: the -1 rotor has no unique logarithm')
	}
	phi := math.atan2(b_norm, s) // half angle, in [0, pi]
	return mv_bivector(phi / b_norm * b[0], phi / b_norm * b[1], phi / b_norm * b[2])
}

// axis_angle extracts the unit axis and rotation angle (in [0, 2*pi]) of a
// unit rotor.
pub fn (r Multivector) axis_angle() ([3]f64, f64) {
	b := r.log()
	bp := b.bivector_part()
	b_norm := math.sqrt(bp[0] * bp[0] + bp[1] * bp[1] + bp[2] * bp[2])
	if b_norm < 1e-12 {
		return [1.0, 0.0, 0.0]!, 0.0
	}
	// b parts are -(t/2)(n3, -n2, n1), so axis = (-b23, b13, -b12)/|b|.
	axis := [-bp[2] / b_norm, bp[1] / b_norm, -bp[0] / b_norm]!
	return axis, 2.0 * b_norm
}

// interpolate slerps between rotors r1 and r2: R(t) = r1 . exp(t . log(r1~ . r2)).
pub fn interpolate(r1 Multivector, r2 Multivector, t f64) Multivector {
	rel := r1.reverse().gp(r2)
	return r1.gp(rel.log().mul_scalar(t).exp())
}

// --- helpers ----------------------------------------------------------------

// popcount counts the set bits of a small non-negative integer.
fn popcount(x int) int {
	mut n := 0
	mut v := x
	for v > 0 {
		n += v & 1
		v = v >> 1
	}
	return n
}

// blade_name returns the symbol of a basis blade, canonical within its grade.
fn blade_name(i int) string {
	return match i {
		0 { '1' }
		1 { 'e1' }
		2 { 'e2' }
		4 { 'e3' }
		3 { 'e12' }
		5 { 'e13' }
		6 { 'e23' }
		7 { 'I' }
		else { '?' }
	}
}
