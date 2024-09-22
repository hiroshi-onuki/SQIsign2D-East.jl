# Julia implementation of SQIsign2D-East

A proof-of-concept implementation of
SQIsign2D-East proposed in
[SQIsign2D-East: A New Signature Scheme Using 2-dimensional Isogenies](https://eprint.iacr.org/2024/771)
using [Julia](https://julialang.org)
and [Nemo](https://www.nemocas.org).
The current implementation is based on the new version
proposed in [Breaking and Repairing SQIsign2D-East](https://eprint.iacr.org/2024/1453),
making it resistant to the attack in that paper.

## Usage

You can run benchmark for SQIsign2D-East and CompactSQIsign2D-East at the NIST security levels 1, 3, and 5 by
```
$ julia bench.jl
```

