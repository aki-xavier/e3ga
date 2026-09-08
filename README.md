# pga

3D Projective Geometric Algebra P(R\*_{3,0,1}) = Cl(3,0,1) for V: 16-component
multivectors, plane-based (dual) PGA — the algebra of Euclidean point / line /
plane geometry, successor of the 8-component Euclidean `e3ga` module.

A multivector has 16 components; every basis blade is indexed by its
generator bitmask (bit 0 = e1, bit 1 = e2, bit 2 = e3, bit 3 = e0):

```
[0]           grade 0:  1
[1][2][4][8]  grade 1:  e1, e2, e3, e0        (planes; e0 = plane at infinity)
[3][5][6]     grade 2:  e12, e13, e23          (Euclidean lines, square -1)
[9][10][12]   grade 2:  e1e0, e2e0, e3e0       (ideal lines, square 0, nilpotent)
[7][11][13][14] grade 3: e123, e12e0, e13e0, e23e0   (points; e123 = origin)
[15]          grade 4:  I = e123e0             (pseudoscalar; I^2 = 0, nilpotent)
```

The metric is degenerate: e1^2 = e2^2 = e3^2 = +1, **e0^2 = 0** — this is what
makes planes, lines and points all live in the same space and makes the outer
product the intersection operator.

## Geometric semantics

| Object | Blade | Representation |
| --- | --- | --- |
| plane | grade 1 | `a e1 + b e2 + c e3 + d e0` for `a x + b y + c z + d = 0` |
| line | grade 2 | direction (e23/e13/e12) + moment (e01/e02/e03) |
| point | grade 3 | `e123 + x e032 + y e013 + z e021` (with e0∧e3∧e2 etc.) |
| motor | even | `M = T . R` — rotation + translation in a single versor |

Native operations:

- **meet (intersection) = outer product**: `plane ^ plane = line`,
  `plane ^ line = point`, point on plane `<=> point ^ plane = 0`.
- **join (union) = regressive product**: `(A* ^ B*)*`, e.g.
  `point ^ ... point v point = line`, three points = plane.
- incidence / distances: `point ^ plane = (n . x + d) . I`.
- motion: `M X M~` applies a motor to any object (points, lines, planes).

## Example

```v
import pga

z := pga.plane([0.0, 0.0, 1.0]!, 0.0)          // the plane z = 0
p := pga.point(1.0, 2.0, 3.0)                    // a point
println(pga.point_plane_dist(p, z))              // 3.0 (signed distance)

l := pga.line_from_points(pga.point(0.0, 0.0, 0.0), pga.point(1.0, 0.0, 0.0))
println(pga.line_angle(l, pga.line_from_points(pga.point(0.0, 0.0, 0.0), pga.point(0.0, 1.0, 0.0)))) // pi/2

m := pga.motor([0.0, 0.0, 1.0]!, math.pi / 2.0, [1.0, 2.0, 3.0]!)
println(m.apply(p).coords())                     // rotate + translate the point
```

## API overview

| Area | Functions |
| --- | --- |
| constructors | mv_zero, mv_scalar, mv_vector (plane coeffs), mv_bivector, mv_trivector, pseudoscalar, e123, e1, e2, e3, e0 |
| access | grade, scalar_part, vector_part, bivector_part, pseudoscalar_part, is_zero, vmax, str |
| arithmetic | add, sub, mul_scalar, div_scalar, neg, eq, copy |
| products | gp, op, ip, lc, rc |
| involutions | reverse, grade_involution, conjugate, dual, undual |
| blade ops | meet (outer product), join (regressive), inverse, norm, normalized |
| motion | rotor, translator, motor, motor_identity, apply, exp, log, interpolate, to_matrix |
| geometry | point, plane, line, line_from_points, coords, ideal_direction, point_dist, point_plane_dist, point_line_dist, line_line_dist, plane_angle, line_angle, project_point_onto_plane, project_point_onto_line, reflect |

## Conventions vs the community

The basis order here is e1-first (bit0 = e1), while bivector.net / ganja.js
use e0-first (`Algebra(3,0,1)`); the component values of each blade differ by
a permutation and sign. The PGA Hodge dual is a blade-wise complement table
(not multiplication by the nilpotent I), converted from PGA4CS Table 4 /
bivector.net to this basis order: `dual(dual(x)) = (-1)^grade(x) . x`,
`undual(x) = (-1)^grade(x) . dual(x)`.

## Testing

```
v test .
```

## History

This module replaces the former 8-component Euclidean GA module (`e3ga`,
Cl(3,0)), which remains in the git history at the commit before the PGA
upgrade.
