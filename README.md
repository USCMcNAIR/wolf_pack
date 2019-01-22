# WOLF PACK

This repository contains **WOLF**ram Language **PACK**ages for the numerical optimization of 1D and 2D structures. In particular, topology optimization of trusses and plane-stress plates.

The aim of WOLF PACK is educational. By performing structural optimization within the interactive environment of a Mathematica notebook, the student can develop computational thinking skills that arise in the design of structures with multiple load-bearing members.

**START DATE**: January-16-2018

**AUTHORS**: `Luis Bahamonde` (luis.bahamonde.jacome@gmail.com), Zafer Gurdal


---
## Summary

The truss topology optimization package (`top1d.wl`) implements, with small variations, the 88 line code proposed by Sokol [1]. It is a standalone package.

The plate topology optimization package (`top2d.wl`) requires the finite element analysis package (`fem2d.wl`), which is also provided in WOLF PACK. Both the `fem2d.wl` and `top2d.wl` packages use [associations](https://reference.wolfram.com/language/ref/Association.html), introduced in version 10.0 of the Wolfram Language.


---
## References

1. T. Sokol, "A 99 line code for discretized Michell truss optimization written in Mathematica," Struct. Multidiscip. Optim., vol. 43, no. 2, pp. 181-190, Feb. 2011.