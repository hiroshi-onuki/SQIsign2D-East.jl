# Julia implementation of SQIsign2D-East

A proof-of-concept implementation of
SQIsign2D-East proposed in
[SQIsign2D-East: A New Signature Scheme Using 2-dimensional Isogenies](https://eprint.iacr.org/2024/****)
using [Julia](https://julialang.org)
and [Nemo](https://www.nemocas.org).

## Usage

You can run benchmark for SQIsign2D-East and CompactSQIsign2D-East at the NIST security levels 1, 3, and 5 by
```
$ julia bench.jl
```

