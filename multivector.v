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

// meet returns the largest common subblade (intersection) of the blades self
// and o.  When their union spans the full 3D space the classic dual formula
// (self* ^ o*)* applies; containment cases return the contained blade, e.g.
// meet(e1, e12) = e1 and meet(e12, e13) = e1.  Non-blade inputs panic.
pub fn (m Multivector) meet(o Multivector) Multivector {
	if m.is_zero() || o.is_zero() {
		return mv_zero()
	}
	ga := m.blade_grade()
	gb := o.blade_grade()
	if ga < 0 || gb < 0 {
		panic('meet: arguments must be blades')
	}
	if ga == 0 || gb == 0 {
		return mv_zero()
	}
	if !m.op(o).is_zero() {
		return mv_zero()
	}
	// The blades share a non-trivial subblade.
	if ga < gb {
		return m
	}
	if gb < ga {
		return o
	}
	// Equal grades: parallel vectors / coincident planes / the pseudoscalar
	// itself share their whole span; distinct planes (possible for grade 2
	// in 3D) meet in their intersection line.
	d := m.dual().op(o.dual())
	if d.is_zero() {
		return m
	}
	return d.dual()
}

// join returns the smallest blade (union span) containing both blades self
// and o: their wedge when they are independent, the larger containing blade
// when they share a subblade, or the pseudoscalar for two distinct planes.
pub fn (m Multivector) join(o Multivector) Multivector {
	if m.is_zero() {
		return o
	}
	if o.is_zero() {
		return m
	}
	ga := m.blade_grade()
	gb := o.blade_grade()
	if ga < 0 || gb < 0 {
		panic('join: arguments must be blades')
	}
	if ga == 0 {
		return o
	}
	if gb == 0 {
		return m
	}
	w := m.op(o)
	if !w.is_zero() {
		return w
	}
	if ga < gb {
		return o
	}
	if gb < ga {
		return m
	}
	if m.dual().op(o.dual()).is_zero() {
		return m
	}
	return pseudoscalar()
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

// inverse returns A^-1 = rev(A)/(A rev(A))_0 for an invertible blade or
// versor; a general multivector with a non-scalar A*rev(A) panics.
pub fn (m Multivector) inverse() Multivector {
	prod := m.gp(m.reverse())
	s := prod.values[0]
	if math.abs(s) < 1e-12 || !prod.grade(0).eq(prod) {
		panic('inverse: only blades and versors with nonzero norm are supported')
	}
	return m.reverse().div_scalar(s)
}

// lc is the left contraction A _| B: sums of <_A_g _B_h>_(h-g) for g <= h.
pub fn (m Multivector) lc(o Multivector) Multivector {
	mut res := Multivector{}
	for ga in 1 .. num_grades {
		a_g := m.grade(ga)
		if a_g.is_zero() {
			continue
		}
		for gb in ga .. num_grades {
			b_g := o.grade(gb)
			if b_g.is_zero() {
				continue
			}
			res = res.add(a_g.gp(b_g).grade(gb - ga))
		}
	}
	return res
}

// rc is the right contraction A |_ B: sums of <_A_g _B_h>_(g-h) for g >= h.
pub fn (m Multivector) rc(o Multivector) Multivector {
	mut res := Multivector{}
	for ga in 1 .. num_grades {
		a_g := m.grade(ga)
		if a_g.is_zero() {
			continue
		}
		for gb in 1 .. ga + 1 {
			b_g := o.grade(gb)
			if b_g.is_zero() {
				continue
			}
			res = res.add(a_g.gp(b_g).grade(ga - gb))
		}
	}
	return res
}

// commutator returns [self, o] = (self o - o self) / 2, the Lie bracket of
// the even subalgebra when both operands are even.
pub fn (m Multivector) commutator(o Multivector) Multivector {
	return m.gp(o).sub(o.gp(m)).mul_scalar(0.5)
}

// anticommutator returns {self, o} = (self o + o self) / 2.
pub fn (m Multivector) anticommutator(o Multivector) Multivector {
	return m.gp(o).add(o.gp(m)).mul_scalar(0.5)
}

// proj projects self onto the blade o: (self . o) o^-1.
pub fn (m Multivector) proj(o Multivector) Multivector {
	return m.ip(o).gp(o.inverse())
}

// rej returns the part of self orthogonal to the blade o.
pub fn (m Multivector) rej(o Multivector) Multivector {
	return m.sub(m.proj(o))
}

// reflect mirrors self across the plane with the given normal:
// v' = -n v n, where n is the normalized normal; works for vectors and blades.
pub fn (m Multivector) reflect(normal [3]f64) Multivector {
	len2 := normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]
	if len2 < 1e-18 {
		panic('reflect: zero normal')
	}
	n := mv_vector(normal[0], normal[1], normal[2]).div_scalar(math.sqrt(len2))
	return n.gp(m).gp(n).neg()
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

// blade_grade returns the grade of a pure blade, or -1 when the multivector
// holds mixed grades or is zero.
fn (m Multivector) blade_grade() int {
	mut g := -1
	for i in 0 .. num_components {
		if m.values[i] != 0.0 {
			gi := popcount(i)
			if g != -1 && gi != g {
				return -1
			}
			g = gi
		}
	}
	return g
}

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
