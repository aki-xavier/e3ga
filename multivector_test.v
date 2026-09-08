module e3ga

// Tests for the plain 3D Euclidean GA Cl(3,0) module.
//
// The geometric product is cross-checked against an independent reference
// implementation (ref_blade_gp) that folds blades back into generator
// vectors and bubble-sorts them, so a systematic sign error in gp_blade
// cannot hide behind a duplicated derivation.
import math

fn almost(a f64, b f64) bool {
	return math.abs(a - b) < 1e-9
}

fn mv_almost(a Multivector, b Multivector) bool {
	for i in 0 .. 8 {
		if math.abs(a.values[i] - b.values[i]) > 1e-9 {
			return false
		}
	}
	return true
}

fn mv_equal(a Multivector, b Multivector) bool {
	for i in 0 .. 8 {
		if math.abs(a.values[i] - b.values[i]) > 1e-6 {
			return false
		}
	}
	return true
}

// ref_blade_gp is the reference basis-blade product: expand both blades into
// generator vectors, bubble-sort the concatenation (each inversion swaps a
// generator pair), and XOR the generator bits.
fn ref_blade_gp(ma int, mb int) (int, f64) {
	mut gens := []int{}
	for i in 0 .. 3 {
		if ma & (1 << i) != 0 {
			gens << i
		}
	}
	for i in 0 .. 3 {
		if mb & (1 << i) != 0 {
			gens << i
		}
	}
	mut sign := 1.0
	mut n := gens.len
	for i in 0 .. n - 1 {
		for j in 0 .. n - 1 - i {
			if gens[j] > gens[j + 1] {
				tmp := gens[j]
				gens[j] = gens[j + 1]
				gens[j + 1] = tmp
				sign = -sign
			}
		}
	}
	mut mask := 0
	for g in gens {
		mask = mask ^ (1 << g)
	}
	return mask, sign
}

fn blade(mask int) Multivector {
	mut m := Multivector{}
	m.values[mask] = 1.0
	return m
}

fn test_basis_squares() {
	assert mv_equal(e1().gp(e1()), mv_scalar(1.0))
	assert mv_equal(e2().gp(e2()), mv_scalar(1.0))
	assert mv_equal(e3().gp(e3()), mv_scalar(1.0))
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).gp(mv_bivector(1.0, 0.0, 0.0)), mv_scalar(-1.0))
	assert mv_equal(mv_bivector(0.0, 1.0, 0.0).gp(mv_bivector(0.0, 1.0, 0.0)), mv_scalar(-1.0))
	assert mv_equal(mv_bivector(0.0, 0.0, 1.0).gp(mv_bivector(0.0, 0.0, 1.0)), mv_scalar(-1.0))
	assert mv_equal(pseudoscalar().gp(pseudoscalar()), mv_scalar(-1.0))
}

fn test_basis_products() {
	assert mv_equal(e1().gp(e2()), mv_bivector(1.0, 0.0, 0.0))
	assert mv_equal(e2().gp(e1()), mv_bivector(-1.0, 0.0, 0.0))
	assert mv_equal(e1().gp(e3()), mv_bivector(0.0, 1.0, 0.0))
	assert mv_equal(e3().gp(e1()), mv_bivector(0.0, -1.0, 0.0))
	assert mv_equal(e2().gp(e3()), mv_bivector(0.0, 0.0, 1.0))
	assert mv_equal(e3().gp(e2()), mv_bivector(0.0, 0.0, -1.0))
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).gp(e3()), pseudoscalar())
	assert mv_equal(e1().gp(mv_bivector(0.0, 0.0, 1.0)), pseudoscalar())
	assert mv_equal(pseudoscalar().gp(e1()), mv_bivector(0.0, 0.0, 1.0))
	assert mv_equal(pseudoscalar().gp(e2()), mv_bivector(0.0, -1.0, 0.0))
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).gp(mv_bivector(0.0, 0.0, 1.0)), mv_bivector(0.0, 1.0, 0.0))
	assert mv_equal(mv_bivector(0.0, 1.0, 0.0).gp(mv_bivector(0.0, 0.0, 1.0)), mv_bivector(-1.0, 0.0, 0.0))
	assert mv_equal(mv_bivector(0.0, 0.0, 1.0).gp(mv_bivector(1.0, 0.0, 0.0)), mv_bivector(0.0, -1.0, 0.0))
}

fn test_gp_reference_crosscheck() {
	for ma in 0 .. 8 {
		for mb in 0 .. 8 {
			dst, sign := gp_blade(ma, mb)
			rdst, rsign := ref_blade_gp(ma, mb)
			assert dst == rdst
			assert almost(sign, rsign)
		}
	}
}

fn test_gp_distributivity_associativity() {
	a := mv_scalar(2.0).add(e1()).add(mv_bivector(0.5, -1.0, 0.25))
	b := mv_vector(1.0, -2.0, 0.5).add(mv_bivector(0.0, 0.25, -0.5))
	c := e3().add(pseudoscalar()).add(mv_scalar(-1.5))
	assert mv_almost(a.add(b).gp(c), a.gp(c).add(b.gp(c)))
	assert mv_almost(a.gp(b).gp(c), a.gp(b.gp(c)))
}

fn test_outer_product() {
	assert mv_equal(e1().op(e1()), mv_zero())
	assert mv_equal(e1().op(e2()), mv_bivector(1.0, 0.0, 0.0))
	assert mv_equal(e2().op(e1()), mv_bivector(-1.0, 0.0, 0.0))
	assert mv_equal(e1().op(e3()), mv_bivector(0.0, 1.0, 0.0))
	assert mv_equal(e2().op(e3()), mv_bivector(0.0, 0.0, 1.0))
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).op(e3()), pseudoscalar())
	assert mv_equal(e1().op(e1().add(e2())), mv_bivector(1.0, 0.0, 0.0))
	// cross-check against the reference for disjoint versus overlapping blades
	for ma in 0 .. 8 {
		for mb in 0 .. 8 {
			mut ref := Multivector{}
			if ma & mb == 0 {
				dst, sign := ref_blade_gp(ma, mb)
				ref.values[dst] = sign
			}
			assert mv_equal(blade(ma).op(blade(mb)), ref)
		}
	}
}

fn test_inner_product() {
	assert mv_equal(e1().ip(e1()), mv_scalar(1.0))
	assert mv_equal(e1().ip(e2()), mv_zero())
	assert mv_equal(e1().ip(mv_scalar(2.0)), mv_zero())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).ip(mv_bivector(1.0, 0.0, 0.0)), mv_scalar(-1.0))
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).ip(mv_bivector(0.0, 0.0, 1.0)), mv_zero())
	assert mv_equal(e1().ip(mv_bivector(1.0, 0.0, 0.0)), e2())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).ip(e1()), e2().neg())
}

fn test_reverse_involution_conjugate() {
	assert mv_equal(mv_scalar(1.0).reverse(), mv_scalar(1.0))
	assert mv_equal(e1().reverse(), e1())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).reverse(), mv_bivector(-1.0, 0.0, 0.0))
	assert mv_equal(pseudoscalar().reverse(), pseudoscalar().neg())
	assert mv_equal(e1().grade_involution(), e1().neg())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).grade_involution(), mv_bivector(1.0, 0.0, 0.0))
	assert mv_equal(pseudoscalar().grade_involution(), pseudoscalar().neg())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).conjugate(), mv_bivector(-1.0, 0.0, 0.0))
}

fn test_dual_undual() {
	assert mv_equal(mv_scalar(1.0).dual(), pseudoscalar().neg())
	assert mv_equal(e1().dual(), mv_bivector(0.0, 0.0, -1.0))
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).dual(), e3())
	assert mv_equal(pseudoscalar().dual(), mv_scalar(1.0))
	assert mv_equal(e1().undual(), mv_bivector(0.0, 0.0, 1.0))
	assert mv_equal(e1().dual().undual(), e1())
}

fn test_meet() {
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(mv_bivector(0.0, 1.0, 0.0)), e1())
}

fn test_norm_normalized() {
	v := mv_vector(1.0, 2.0, 2.0)
	assert almost(v.norm(), 3.0)
	u := v.normalized()
	want := mv_vector(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)
	assert mv_equal(u, want)
	assert almost(mv_zero().normalized().norm(), 0.0)
}

fn test_rotor_axis_z() {
	r := rotor([0.0, 0.0, 1.0]!, math.pi / 2.0)
	assert almost(r.scalar_part(), math.cos(math.pi / 4.0))
	got1 := r.apply(e1())
	want1 := mv_vector(0.0, 1.0, 0.0)
	got2 := r.apply(e2())
	want2 := mv_vector(-1.0, 0.0, 0.0)
	got3 := r.apply(e3())
	assert mv_almost(got1, want1)
	assert mv_almost(got2, want2)
	assert mv_almost(got3, e3())
	// axis by magnitude: rotor of angle theta rotates so that e1 -> cos e1 + sin e2
	r30 := rotor([0.0, 0.0, 1.0]!, math.pi / 6.0)
	w := r30.apply(e1())
	assert almost(w.vector_part()[0], math.cos(math.pi / 6.0))
	assert almost(w.vector_part()[1], math.sin(math.pi / 6.0))
}

fn test_rotor_axis_y() {
	r := rotor([0.0, 1.0, 0.0]!, math.pi / 2.0)
	assert mv_almost(r.apply(e1()), mv_vector(0.0, 0.0, -1.0))
}

fn test_rotor_preserves_norm_dot() {
	r := rotor([1.0, 2.0, 3.0]!, 1.234)
	v := mv_vector(0.4, -1.2, 0.7)
	assert almost(r.apply(v).norm(), v.norm())
	// axis itself is invariant under its own rotation
	axis := mv_vector(1.0, 2.0, 3.0).normalized()
	assert mv_almost(r.apply(axis), axis)
}

fn test_exp_log_roundtrip() {
	r := rotor([1.0, 2.0, 3.0]!, 1.234)
	back := r.log().exp()
	assert mv_almost(r, back)
	b := mv_bivector(0.3, -0.2, 0.5)
	e := b.exp()
	back_b := e.log()
	assert mv_almost(b, back_b)
	// scalar + bivector handled by the closed form: exp(s + B) = e^s . exp(B)
	sb := mv_scalar(0.7).add(b)
	assert mv_almost(sb.exp(), b.exp().mul_scalar(math.exp(0.7)))
}

fn test_interpolate() {
	r1 := rotor([0.0, 0.0, 1.0]!, 0.0)
	r2 := rotor([0.0, 0.0, 1.0]!, math.pi / 2.0)
	mid := interpolate(r1, r2, 0.5)
	w := mid.apply(e1())
	assert almost(w.vector_part()[0], math.cos(math.pi / 4.0))
	assert almost(w.vector_part()[1], math.sin(math.pi / 4.0))
	assert mv_almost(interpolate(r1, r2, 1.0), r2)
	assert mv_almost(interpolate(r1, r2, 0.0), r1)
}

fn test_axis_angle() {
	r := rotor([0.0, 0.0, 1.0]!, 2.0)
	axis, angle := r.axis_angle()
	assert almost(axis[0], 0.0)
	assert almost(axis[1], 0.0)
	assert almost(axis[2], 1.0)
	assert almost(angle, 2.0)
	r2 := rotor([1.0, 0.0, 0.0]!, math.pi)
	axis2, angle2 := r2.axis_angle()
	assert almost(axis2[0], 1.0)
	assert almost(angle2, math.pi)
}

fn test_operators_str() {
	a := mv_vector(1.0, 2.0, 3.0)
	b := mv_vector(4.0, -1.0, 0.0)
	assert mv_equal(a.add(b), mv_vector(5.0, 1.0, 3.0))
	assert mv_equal(a.sub(b), mv_vector(-3.0, 3.0, 3.0))
	assert mv_equal(a.mul_scalar(2.0), mv_vector(2.0, 4.0, 6.0))
	assert mv_equal(a.div_scalar(2.0), mv_vector(0.5, 1.0, 1.5))
	assert mv_equal(a.neg(), mv_vector(-1.0, -2.0, -3.0))
	assert a.eq(a.copy())
	s := e1().add(mv_bivector(1.0, 0.0, 0.0)).str()
	assert s.contains('e1')
	assert s.contains('e12')
	assert mv_zero().str().contains('0')
	assert almost(pseudoscalar().pseudoscalar_part(), 1.0)
}

fn test_join() {
	assert mv_equal(e1().join(e2()), mv_bivector(1.0, 0.0, 0.0))
	// containment: join is the larger blade
	assert mv_equal(e1().join(mv_bivector(1.0, 0.0, 0.0)), mv_bivector(1.0, 0.0, 0.0))
	// independent line and plane span the whole space
	assert mv_equal(e1().join(mv_bivector(0.0, 0.0, 1.0)), pseudoscalar())
	// distinct planes through a common line span the whole space
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).join(mv_bivector(0.0, 1.0, 0.0)), pseudoscalar())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).join(mv_bivector(0.0, 0.0, 1.0)), pseudoscalar())
	// distinct plane that shares only a line with e12 (2e12 + e23)
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).join(mv_bivector(2.0, 0.0, 1.0)), pseudoscalar())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).join(mv_bivector(1.0, 0.0, 0.0)), mv_bivector(1.0, 0.0, 0.0))
	// parallel vectors share the same span
	assert mv_equal(e2().join(e2().mul_scalar(-3.0)), e2())
	// scalar / zero blades join as the empty blade
	assert mv_equal(mv_scalar(2.0).join(e1()), e1())
	assert mv_equal(mv_zero().join(e1()), e1())
	assert mv_equal(pseudoscalar().join(e1()), pseudoscalar())
}

fn test_meet_containment() {
	// containment cases that the dual formula alone gets wrong
	assert mv_equal(e1().meet(mv_bivector(1.0, 0.0, 0.0)), e1())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(e1()), e1())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(mv_bivector(1.0, 0.0, 0.0)), mv_bivector(1.0, 0.0, 0.0))
	assert mv_equal(e1().meet(pseudoscalar()), e1())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(pseudoscalar()), mv_bivector(1.0, 0.0, 0.0))
	// disjoint blades meet in the empty blade
	assert mv_equal(e1().meet(e2()), mv_zero())
	assert mv_equal(e1().meet(mv_bivector(0.0, 0.0, 1.0)), mv_zero())
	assert mv_equal(mv_scalar(3.0).meet(e1()), mv_zero())
	// distinct planes meet in their intersection line
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(mv_bivector(0.0, 1.0, 0.0)), e1())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(mv_bivector(0.0, 0.0, 1.0)), e2())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).meet(mv_bivector(2.0, 0.0, 1.0)), e2())
}

fn test_inverse() {
	// blade inverse: A rev(A) = |A|^2, so e12^-1 = -e12
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).inverse(), mv_bivector(-1.0, 0.0, 0.0))
	v := mv_vector(1.0, 2.0, 2.0)
	assert mv_equal(v.gp(v.inverse()), mv_scalar(1.0))
	r := rotor([1.0, 2.0, 3.0]!, 1.234)
	assert mv_equal(r.gp(r.inverse()), mv_scalar(1.0))
	assert mv_equal(r.inverse(), r.reverse())
}

fn test_projection_rejection() {
	plane := mv_bivector(1.0, 0.0, 0.0)
	v := mv_vector(1.0, 2.0, 3.0)
	assert mv_equal(v.proj(plane), mv_vector(1.0, 2.0, 0.0))
	assert mv_equal(v.rej(plane), mv_vector(0.0, 0.0, 3.0))
	assert mv_equal(v.proj(plane).add(v.rej(plane)), v)
	assert mv_equal(mv_vector(1.0, 2.0, 3.0).proj(e1()), mv_vector(1.0, 0.0, 0.0))
	assert mv_equal(e3().proj(plane), mv_zero())
}

fn test_reflect() {
	v := mv_vector(1.0, 2.0, 3.0)
	assert mv_equal(v.reflect([0.0, 0.0, 1.0]!), mv_vector(1.0, 2.0, -3.0))
	assert mv_equal(v.reflect([1.0, 0.0, 0.0]!), mv_vector(-1.0, 2.0, 3.0))
	// the normal is normalized internally
	assert mv_equal(v.reflect([0.0, 0.0, 2.0]!), mv_vector(1.0, 2.0, -3.0))
	assert mv_equal(e3().reflect([0.0, 0.0, 1.0]!), e3().neg())
}

fn test_commutator_anticommutator() {
	assert mv_equal(e1().commutator(e2()), mv_bivector(1.0, 0.0, 0.0))
	assert mv_equal(e1().commutator(e1()), mv_zero())
	assert mv_equal(e1().anticommutator(e2()), mv_zero())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).commutator(mv_bivector(0.0, 0.0, 1.0)), mv_bivector(0.0, 1.0, 0.0))
}

fn test_contractions() {
	assert mv_equal(e1().lc(mv_bivector(1.0, 0.0, 0.0)), e2())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).lc(e1()), mv_zero())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).rc(e1()), e2().neg())
	assert mv_equal(e1().rc(mv_bivector(1.0, 0.0, 0.0)), mv_zero())
	assert mv_equal(mv_bivector(1.0, 0.0, 0.0).lc(mv_bivector(0.0, 0.0, 1.0)), mv_zero())
}
