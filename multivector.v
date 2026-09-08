module pga

// 3D Projective Geometric Algebra P(R*_{3,0,1}) = Cl(3,0,1) — the
// "plane-based" (dual) PGA: the algebra of Euclidean 3D point/line/plane
// geometry, successor of the 8-component Euclidean module.
//
// Basis generators {e1, e2, e3, e0} with e1^2 = e2^2 = e3^2 = +1, e0^2 = 0,
// generators anticommute pairwise.  Every basis blade is indexed by its
// bitmask (bit 0 = e1, bit 1 = e2, bit 2 = e3, bit 3 = e0), so the geometric
// product of two basis blades is a XOR/sign computation with an extra metric
// factor when both factors contain the null generator e0:
//
//   [0] grade 0: 1                        [8]  grade 1: e0    (plane at infinity)
//   [1][2][4]   grade 1: e1, e2, e3       [9][10][12] grade 2: e1e0, e2e0, e3e0 (ideal lines)
//   [3][5][6]   grade 2: e12, e13, e23    [11][13][14] grade 3: e12e0, e13e0, e23e0
//   [7] grade 3: e123 (the origin point)  [15] grade 4: I = e123e0 (null, I^2 = 0)
//
// Geometric semantics (the reason for the degenerate metric):
//   plane = grade-1 vector  n + d e0         (n unit normal, distance d)
//   line  = grade-2 blade   direction + moment (Euclidean, squares negative)
//   point = grade-3 blade   e123 + x e032 + y e013 + z e021
//   meet  = outer product (plane ^ plane = line, plane ^ line = point)
//   join  = regressive product (point v point = line, three points = plane)
//   motor  = even-grade versor (rotation + translation), M X M~ applies it.
//
// The pseudoscalar I = e123 e0 squares to 0 (nilpotent), so the Hodge dual is
// NOT multiplication by I; it is the blade-wise complement below (PGA4CS
// Table 4 / bivector.net, converted to this basis order), and
// dual(dual(x)) = (-1)^grade(x) . x.
import math

pub const num_components = 16
pub const num_grades = 5

// Metric of the four generators in bit order (e1, e2, e3, e0).
const generator_metric = [1.0, 1.0, 1.0, 0.0]

pub struct Multivector {
pub mut:
	values [16]f64
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

// mv_vector builds a plane (grade-1) from coefficients (a, b, c, d): the
// plane a x + b y + c z + d = 0.
pub fn mv_vector(a f64, b f64, c f64, d f64) Multivector {
	mut m := Multivector{}
	m.values[1] = a
	m.values[2] = b
	m.values[4] = c
	m.values[8] = d
	return m
}

// mv_bivector builds a bivector from the six components in the order
// (e12, e13, e23, e1e0, e2e0, e3e0).
pub fn mv_bivector(b12 f64, b13 f64, b23 f64, b01 f64, b02 f64, b03 f64) Multivector {
	mut m := Multivector{}
	m.values[3] = b12
	m.values[5] = b13
	m.values[6] = b23
	m.values[9] = b01
	m.values[10] = b02
	m.values[12] = b03
	return m
}

// mv_trivector builds a trivector from the four components in the order
// (e123, e12e0, e13e0, e23e0) — useful for raw point-like blades.
pub fn mv_trivector(t123 f64, t120 f64, t130 f64, t230 f64) Multivector {
	mut m := Multivector{}
	m.values[7] = t123
	m.values[11] = t120
	m.values[13] = t130
	m.values[14] = t230
	return m
}

// pseudoscalar returns I = e123 e0 (squares to 0).
pub fn pseudoscalar() Multivector {
	mut m := Multivector{}
	m.values[15] = 1.0
	return m
}

// e123 returns I3 = e1 e2 e3, the Euclidean pseudoscalar (squares to -1).
pub fn e123() Multivector {
	mut m := Multivector{}
	m.values[7] = 1.0
	return m
}

// Basis vectors (grade 1).  Fresh values per call.
pub fn e1() Multivector {
	return mv_vector(1.0, 0.0, 0.0, 0.0)
}

pub fn e2() Multivector {
	return mv_vector(0.0, 1.0, 0.0, 0.0)
}

pub fn e3() Multivector {
	return mv_vector(0.0, 0.0, 1.0, 0.0)
}

// e0 is the plane at infinity.
pub fn e0() Multivector {
	return mv_vector(0.0, 0.0, 0.0, 1.0)
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

// vector_part returns the four grade-1 components (e1, e2, e3, e0).
pub fn (m Multivector) vector_part() [4]f64 {
	return [m.values[1], m.values[2], m.values[4], m.values[8]]!
}

// bivector_part returns the six grade-2 components in the order
// (e12, e13, e23, e1e0, e2e0, e3e0).
pub fn (m Multivector) bivector_part() [6]f64 {
	return [m.values[3], m.values[5], m.values[6], m.values[9], m.values[10], m.values[12]]!
}

// pseudoscalar_part returns the grade-4 component.
pub fn (m Multivector) pseudoscalar_part() f64 {
	return m.values[15]
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
// sign counts inversions over the generator list; a shared e0 factor makes
// the term vanish because e0^2 = 0.
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
			dst, coeff := gp_blade(ma, mb)
			if coeff == 0.0 {
				continue
			}
			res.values[dst] += coeff * a * b
		}
	}
	return res
}

// gp_blade returns the product of two basis blades as (result mask, coeff).
fn gp_blade(ma int, mb int) (int, f64) {
	mut swaps := 0
	for i in 0 .. 4 {
		if mb & (1 << i) != 0 {
			swaps += popcount(ma >> (i + 1))
		}
	}
	// common basis factors square to their metric; the null generator e0
	// squares to 0, so any shared e0 factor annihilates the term.
	mut coeff := 1.0
	for i in 0 .. 4 {
		if (ma & mb) & (1 << i) != 0 {
			coeff *= generator_metric[i]
		}
	}
	mut sign := 1.0
	if swaps % 2 != 0 {
		sign = -1.0
	}
	return ma ^ mb, sign * coeff
}

// op computes the outer product a ^ b (the PGA meet): zero for blades
// sharing a basis factor, else +/- the union blade.
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
			dst, coeff := gp_blade(ma, mb)
			res.values[dst] += coeff * a * b
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

// --- Hodge dual and regressive joins ----------------------------------------

// dual is the PGA Hodge / Poincare complement: blade-wise with the signs of
// PGA4CS Table 4 (bivector.net), converted to the (e1,e2,e3,e0) bit order.
// dual(dual(x)) = (-1)^grade(x) . x.
pub fn (m Multivector) dual() Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		if m.values[i] != 0.0 {
			res.values[dual_dst[i]] += dual_sign[i] * m.values[i]
		}
	}
	return res
}

// undual returns (-1)^grade(x) . dual(x), the inverse of dual.
pub fn (m Multivector) undual() Multivector {
	mut res := Multivector{}
	for i in 0 .. num_components {
		if m.values[i] != 0.0 {
			sgn := dual_sign[i]
			if popcount(dual_dst[i]) % 2 != 0 {
				// (-1)^grade of the result: odd grades flip
				res.values[dual_dst[i]] += -sgn * m.values[i]
			} else {
				res.values[dual_dst[i]] += sgn * m.values[i]
			}
		}
	}
	return res
}

// meet returns the intersection of two blades: the outer product (native PGA
// incidence), e.g. meet(plane, plane) = line, meet(plane, line) = point.
pub fn (m Multivector) meet(o Multivector) Multivector {
	return m.op(o)
}

// join returns the union of two blades via the regressive product
// (self* ^ o*)*, e.g. join(point, point) = line, join(line, point) = plane.
pub fn (m Multivector) join(o Multivector) Multivector {
	return m.dual().op(o.dual()).dual()
}

// --- norms and inverses -----------------------------------------------------

// norm returns sqrt(|<self.reverse(self)>_0|); lines and rotors have norm 1,
// null blades (ideal lines, points) have norm 0.
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

// inverse returns A^-1 = rev(A) (s - p I)/s^2 with A rev(A) = s + p I — exact
// for blades and motors (unit elements give rev(A)); null elements (points,
// ideal lines) panic.
pub fn (m Multivector) inverse() Multivector {
	prod := m.gp(m.reverse())
	s := prod.values[0]
	p := prod.values[15]
	if math.abs(s) < 1e-12 {
		panic('inverse: degenerate null element')
	}
	// (s + p I)^-1 = (s - p I) / s^2 since I^2 = 0
	inv := m.reverse().mul_scalar(s).sub(m.reverse().gp(pseudoscalar()).mul_scalar(p))
	return inv.div_scalar(s * s)
}

// --- motion (the even subalgebra: dual quaternions of SE(3)) -----------------

// rotor returns the rotation rotor by `angle` radians about the axis through
// the origin in direction `axis`: R = exp(-angle/2 . L) with L = unit line.
pub fn rotor(axis [3]f64, angle f64) Multivector {
	len2 := axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]
	if len2 < 1e-18 {
		panic('rotor: zero axis')
	}
	inv_len := 1.0 / math.sqrt(len2)
	u := mv_vector(axis[0] * inv_len, axis[1] * inv_len, axis[2] * inv_len, 0.0)
	// unit line through the origin: u _| I3
	line := u.lc(e123())
	half := angle / 2.0
	return mv_scalar(math.cos(half)).sub(line.mul_scalar(math.sin(half)))
}

// translator returns the translator T = exp(-d/2 . e0 ^ t) = 1 - (e0 ^ t)/2
// (the ideal line e0 ^ t is nilpotent, so the series truncates).
pub fn translator(displacement [3]f64) Multivector {
	t := mv_vector(displacement[0], displacement[1], displacement[2], 0.0)
	return mv_scalar(1.0).sub(e0().op(t).mul_scalar(0.5))
}

// motor builds the rigid motion M = T . R (rotate then translate).
pub fn motor(axis [3]f64, angle f64, displacement [3]f64) Multivector {
	return translator(displacement).gp(rotor(axis, angle))
}

// motor_identity returns the identity motor.
pub fn motor_identity() Multivector {
	return mv_scalar(1.0)
}

// apply returns M . v . M~ for a versor M; rotates/translates points, lines
// and planes.
pub fn (m Multivector) apply(v Multivector) Multivector {
	return m.gp(v).gp(m.reverse())
}

// exp exponentiates a scalar + bivector to the motor group: exp(s + B) =
// e^s exp(B).  Simple Euclidean bivectors use cos/sin, ideal (nilpotent)
// bivectors truncate at 1 + B, and general screw bivectors use the axis-pair
// closed form.
pub fn (m Multivector) exp() Multivector {
	if !m.grade(1).is_zero() || !m.grade(3).is_zero() || !m.grade(4).is_zero() {
		panic('exp: only scalar + bivector are supported')
	}
	scale := math.exp(m.values[0])
	b := m.grade(2)
	if b.is_zero() {
		return mv_scalar(scale)
	}
	b2 := b.gp(b)
	s := b2.values[0]
	if math.abs(s) < 1e-12 {
		// nilpotent: B^2 = 0 (pure translation / ideal line)
		return mv_scalar(scale).add(b.mul_scalar(scale))
	}
	p := b2.values[15]
	u := math.sqrt(-s)
	v := -p / (2.0 * u)
	perp := b.gp(pseudoscalar()).div_scalar(u) // theta_hat _| = B.I / u
	hat := b.sub(perp.mul_scalar(v)).div_scalar(u)
	cu := scale * math.cos(u)
	su := scale * math.sin(u)
	return mv_scalar(cu).add(hat.mul_scalar(su)).gp(mv_scalar(1.0).add(perp.mul_scalar(v)))
}

// log returns the bivector B with exp(B) = self for a unit motor.  Pure
// translations return their ideal line; the -1 motor panics.
pub fn (m Multivector) log() Multivector {
	a := m.values[0]
	b := m.grade(2)
	c := m.values[15]
	sin_u := math.sqrt(math.max(0.0, 1.0 - a * a))
	if sin_u < 1e-9 {
		if b.is_zero() {
			if a < 0.0 {
				panic('log: the -1 motor has no unique logarithm')
			}
			return mv_zero()
		}
		return b // pure translation
	}
	u := math.atan2(sin_u, a)
	v := -c / sin_u
	perp := b.gp(pseudoscalar()).div_scalar(sin_u) // sin_u . theta_hat _| = B.I
	beta := v * math.cos(u)
	hat := b.sub(perp.mul_scalar(beta)).div_scalar(sin_u)
	return hat.mul_scalar(u).add(perp.mul_scalar(v))
}

// interpolate slerps between motors m1 and m2: M(t) = m1 . exp(t . log(m1~ . m2)).
pub fn interpolate(m1 Multivector, m2 Multivector, t f64) Multivector {
	rel := m1.reverse().gp(m2)
	return m1.gp(rel.log().mul_scalar(t).exp())
}

// to_matrix returns the equivalent 4x4 homogeneous transform [R|t], flattened
// row-major into 16 components (row r, col c is at index 4*r + c).
pub fn (m Multivector) to_matrix() [16]f64 {
	origin_t := m.apply(point(0.0, 0.0, 0.0))
	px := m.apply(point(1.0, 0.0, 0.0)).coords()
	py := m.apply(point(0.0, 1.0, 0.0)).coords()
	pz := m.apply(point(0.0, 0.0, 1.0)).coords()
	tx := origin_t.coords()
	return [
		px[0] - tx[0],
		py[0] - tx[0],
		pz[0] - tx[0],
		tx[0],
		px[1] - tx[1],
		py[1] - tx[1],
		pz[1] - tx[1],
		tx[1],
		px[2] - tx[2],
		py[2] - tx[2],
		pz[2] - tx[2],
		tx[2],
		0.0,
		0.0,
		0.0,
		1.0,
	]!
}

// --- blade helpers ----------------------------------------------------------

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

// blade_name returns the symbol of a basis blade, sorted by generator
// (e1, e2, e3, e0).
fn blade_name(i int) string {
	return match i {
		0 { '1' }
		1 { 'e1' }
		2 { 'e2' }
		4 { 'e3' }
		8 { 'e0' }
		3 { 'e12' }
		5 { 'e13' }
		6 { 'e23' }
		9 { 'e1e0' }
		10 { 'e2e0' }
		12 { 'e3e0' }
		7 { 'e123' }
		11 { 'e12e0' }
		13 { 'e13e0' }
		14 { 'e23e0' }
		15 { 'I' }
		else { '?' }
	}
}

// The Hodge complement: blade -> (complement slot, sign).  Converted from
// PGA4CS Table 4 / bivector.net: 1<->I, e0<->e123, e1<->e032, e2<->e013,
// e3<->e021, e01<->e23, e02<->e31, e03<->e12, e032->-e1, e013->-e2,
// e021->-e3, e123->-e0.
const dual_dst = [
	15, // 1    -> I
	14, // e1   -> e032 (= e23e0 slot)
	13, // e2   -> e013 (= e13e0 slot)
	12, // e12  -> e03 (= e3e0 slot)
	11, // e3   -> e021 (= e12e0 slot)
	10, // e13  -> -e02 (= e2e0 slot)
	9, // e23  -> e01 (= e1e0 slot)
	8, // e123 -> -e0
	7, // e0   -> e123
	6, // e01  -> e23
	5, // e02  -> -e31 (= e13 slot, negated)
	4, // e021 -> -e3
	3, // e03  -> e12
	2, // e013 -> -e2
	1, // e032 -> -e1
	0, // I    -> 1
]

const dual_sign = [
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	-1.0,
	1.0,
	-1.0,
	1.0,
	1.0,
	-1.0,
	-1.0,
	1.0,
	-1.0,
	-1.0,
	1.0,
]
