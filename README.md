# e3ga

Plain 3D Euclidean geometric algebra Cl(3,0) for V: 8-component
multivectors over f64, blades, rotors, and the geometric object operations
built on them.

A multivector has 8 components; each basis blade is indexed by its
basis-vector mask (bit 0 = e1, bit 1 = e2, bit 2 = e3):

```
[0] grade 0:  1
[1] [2] [4]   grade 1: e1, e2, e3
[3] [5] [6]   grade 2: e12, e13, e23
[7]           grade 3: I = e123
```

The metric is Euclidean (e1^2 = e2^2 = e3^2 = +1, ei ej = -ej ei for i != j).
The even subalgebra (scalar + bivector) is isomorphic to the quaternions:
a unit rotor R = cos(t/2) - sin(t/2) B-hat rotates by v' = R v R-twiddle.

## Features

- Geometric, outer, and Hestenes fat-dot products; left/right contractions;
  commutator and anticommutator
- Reversal, grade involution, Clifford conjugate, Hodge dual (dual/undual)
- meet and join on blades (containment cases included), projection/rejection,
  reflection across a plane
- Blades and versors inverses
- Rotors: axis-angle construction, sandwich application (apply), exp/log of
  even multivectors, axis-angle extraction, slerp interpolation

## Quick example

```v
import e3ga
import math

v := e3ga.mv_vector(1.0, 2.0, 3.0)
plane := e3ga.mv_bivector(1.0, 0.0, 0.0) // the xy plane
println(v.proj(plane).str())             // Multivector(+1.0000*e1 +2.0000*e2)

r := e3ga.rotor([0.0, 0.0, 1.0]!, math.pi / 2.0)
println(r.apply(e3ga.e1()).str())        // Multivector(+1.0000*e2)

w := e3ga.mv_vector(0.0, 1.0, 0.0)
println(w.meet(plane).str())             // Multivector(+1.0000*e2)

println(v.reflect([0.0, 0.0, 1.0]!).str()) // Multivector(+1.0000*e1 +2.0000*e2 -3.0000*e3)
```

## API overview

| Area | Functions |
| --- | --- |
| constructors | mv_zero, mv_scalar, mv_vector, mv_bivector, pseudoscalar, e1, e2, e3 |
| access | grade, scalar_part, vector_part, bivector_part, pseudoscalar_part, is_zero, vmax, str |
| arithmetic | add, sub, mul_scalar, div_scalar, neg, eq, copy |
| products | gp, op, ip, lc, rc, commutator, anticommutator |
| involutions | reverse, grade_involution, conjugate, dual, undual |
| blade ops | meet, join, proj, rej, inverse, reflect, blade_grade |
| norms | norm, normalized |
| rotors | rotor, apply, exp, log, axis_angle, interpolate |

`join` and `meet` operate on blades and panic on non-blade inputs;
`inverse` covers blades and versors (elements with a scalar A*rev(A)).
`exp`/`log` handle even multivectors (scalar + bivector) only.

## Testing

```
v test .
```
