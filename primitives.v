module pga

// Geometry layer for PGA: points, planes, lines and the distance / angle /
// projection toolbox.
//
//   plane:  n + d e0            (unit normal n, signed distance d)
//   line:   bivector direction + moment (join of two points; Euclidean lines
//           square to -|d|^2, ideal lines are null)
//   point:  trivector e123 + x e032 + y e013 + z e021
//
// Native incidence: point on plane  <=> point ^ plane = 0
//                   point on line   <=> line ^ plane-at-? (line ^ point = 0)
//                   plane ^ plane   = line;  plane ^ line = point.
import math

// point returns the PGA point at (x, y, z).
pub fn point(x f64, y f64, z f64) Multivector {
	mut p := e123()
	p = p.add(e0().op(e3()).op(e2()).mul_scalar(x))
	p = p.add(e0().op(e1()).op(e3()).mul_scalar(y))
	p = p.add(e0().op(e2()).op(e1()).mul_scalar(z))
	return p
}

// plane returns the unit plane n . x + d = 0 (normal normalized internally).
pub fn plane(normal [3]f64, d f64) Multivector {
	nl := math.sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2])
	if nl < 1e-12 {
		panic('plane normal is zero')
	}
	return mv_vector(normal[0] / nl, normal[1] / nl, normal[2] / nl, d)
}

// line_from_points returns the join of two points: the line through them
// (the ideal line if the points coincide in direction).
pub fn line_from_points(p1 Multivector, p2 Multivector) Multivector {
	return p1.join(p2)
}

// line returns the line through `p` in direction `dir` (direction need not be
// unit; the vanish point of `dir` is joined with p).
pub fn line(p Multivector, dir [3]f64) Multivector {
	u := mv_vector(dir[0], dir[1], dir[2], 0.0)
	vanishing := e0().op(u.lc(e123()))
	return p.join(vanishing)
}

// coords extracts the euclidean coordinates (x, y, z) of a trivector point
// (normalising by the e123 weight, so scaled or reflected points work).
pub fn (m Multivector) coords() [3]f64 {
	w := m.values[7]
	if math.abs(w) < 1e-9 {
		panic('multivector has no e123 component; not a finite point')
	}
	return [-m.values[14] / w, m.values[13] / w, -m.values[11] / w]!
}

// ideal_direction returns the direction of a line: the direction part of its
// vanishing point (line ^ e0), which is the ideal point of the line.
pub fn (m Multivector) ideal_direction() [3]f64 {
	v := m.meet(e0())
	return [-v.values[14], v.values[13], -v.values[11]]!
}

// point_on_line returns either a finite point of `line` (intersecting it with
// the first axis plane that is not parallel to it) or the ideal point.
pub fn (line Multivector) point_on_line() Multivector {
	for normal in [[1.0, 0.0, 0.0]!, [0.0, 1.0, 0.0]!, [0.0, 0.0, 1.0]!] {
		q := line.meet(plane(normal, 0.0))
		if !q.is_zero() && math.abs(q.values[7]) > 1e-6 {
			return q
		}
	}
	return line.meet(e0())
}

// --- distances ---------------------------------------------------------------

// point_dist returns the euclidean distance between two points.
pub fn point_dist(a Multivector, b Multivector) f64 {
	c1 := a.coords()
	c2 := b.coords()
	dx := c1[0] - c2[0]
	dy := c1[1] - c2[1]
	dz := c1[2] - c2[2]
	return math.sqrt(dx * dx + dy * dy + dz * dz)
}

// point_plane_dist returns the signed distance of point p to plane pi
// (positive on the side of the normal).
pub fn point_plane_dist(p Multivector, pi Multivector) f64 {
	return p.meet(pi).values[15]
}

// point_line_dist returns the distance from point p to the line l.
pub fn point_line_dist(p Multivector, l Multivector) f64 {
	q := l.point_on_line()
	d := l.ideal_direction()
	an := math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
	if an < 1e-12 {
		return point_dist(p, q)
	}
	c := p.coords()
	// vector from q to p
	rx := c[0] - q.coords()[0]
	ry := c[1] - q.coords()[1]
	rz := c[2] - q.coords()[2]
	// distance = |r x d| / |d|
	crx := ry * d[2] - rz * d[1]
	cry := rz * d[0] - rx * d[2]
	crz := rx * d[1] - ry * d[0]
	return math.sqrt(crx * crx + cry * cry + crz * crz) / an
}

// line_line_dist returns the distance between two lines (0 for intersecting
// lines).
pub fn line_line_dist(a Multivector, b Multivector) f64 {
	da := a.ideal_direction()
	db := b.ideal_direction()
	na := math.sqrt(da[0] * da[0] + da[1] * da[1] + da[2] * da[2])
	nb := math.sqrt(db[0] * db[0] + db[1] * db[1] + db[2] * db[2])
	if na < 1e-12 || nb < 1e-12 {
		return 1e300
	}
	// common normal direction
	nx := da[1] * db[2] - da[2] * db[1]
	ny := da[2] * db[0] - da[0] * db[2]
	nz := da[0] * db[1] - da[1] * db[0]
	nn := math.sqrt(nx * nx + ny * ny + nz * nz)
	if nn < 1e-12 {
		// parallel lines: distance from a point of b to line a
		return point_line_dist(b.point_on_line(), a)
	}
	pa := a.point_on_line()
	pb := b.point_on_line()
	ca := pa.coords()
	cb := pb.coords()
	// (b - a) . n / |n|
	dx := cb[0] - ca[0]
	dy := cb[1] - ca[1]
	dz := cb[2] - ca[2]
	return math.abs(dx * nx + dy * ny + dz * nz) / nn
}

// --- angles ------------------------------------------------------------------

// plane_angle returns the angle (radians, [0, pi/2]) between two planes.
pub fn plane_angle(a Multivector, b Multivector) f64 {
	mut dot := math.abs(a.ip(b).values[0])
	if dot > 1.0 {
		dot = 1.0
	}
	return math.acos(dot)
}

// line_angle returns the angle (radians, [0, pi/2]) between two lines.
pub fn line_angle(a Multivector, b Multivector) f64 {
	da := a.ideal_direction()
	db := b.ideal_direction()
	na := math.sqrt(da[0] * da[0] + da[1] * da[1] + da[2] * da[2])
	nb := math.sqrt(db[0] * db[0] + db[1] * db[1] + db[2] * db[2])
	if na < 1e-12 || nb < 1e-12 {
		return 0.0
	}
	mut dot := (da[0] * db[0] + da[1] * db[1] + da[2] * db[2]) / (na * nb)
	if dot < -1.0 {
		dot = -1.0
	}
	if dot > 1.0 {
		dot = 1.0
	}
	return math.acos(math.abs(dot))
}

// --- projections -------------------------------------------------------------

// project_point_onto_plane projects p onto the plane pi along its normal.
pub fn project_point_onto_plane(p Multivector, pi Multivector) Multivector {
	pn := pi.vector_part()
	c := p.coords()
	dist := pn[0] * c[0] + pn[1] * c[1] + pn[2] * c[2] + pi.values[8]
	return point(c[0] - dist * pn[0], c[1] - dist * pn[1], c[2] - dist * pn[2])
}

// project_point_onto_line projects p onto the line l.
pub fn project_point_onto_line(p Multivector, l Multivector) Multivector {
	q := l.point_on_line()
	d := l.ideal_direction()
	an := math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
	if an < 1e-12 {
		return q
	}
	c := p.coords()
	qc := q.coords()
	lam := ((c[0] - qc[0]) * d[0] + (c[1] - qc[1]) * d[1] + (c[2] - qc[2]) * d[2]) / (an * an)
	return point(qc[0] + lam * d[0], qc[1] + lam * d[1], qc[2] + lam * d[2])
}

// reflect mirrors self across the plane with the given normal: X' = p X p
// with p the unit plane vector.
pub fn (m Multivector) reflect(normal [3]f64) Multivector {
	p := plane(normal, 0.0)
	return p.gp(m).gp(p)
}
