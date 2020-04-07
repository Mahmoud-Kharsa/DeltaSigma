# DeltaSigma

A toolbox for delta-sigma modulator design in Julia. Adapted from the MATLAB
[Delta-Sigma Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox)

## Installation

To install enter the Pkg REPL by pressing `]` from the Julia REPL and enter:

```julia
add https://github.com/Mahmoud-Kharsa/DeltaSigma.git
```

## Usage

The toolbox offers functions similar to the ones available in the MATLAB toolbox.
All functions have docstrings, which can be viewed from the REPL, for example
`?synthesizeNTF`

### Example

Here we simulate delta sigma modulation with a 3<sup>rd</sup> order modulator
on a sine wave input.

```julia
using DeltaSigma
using PyPlot

H = synthesizeNTF(3, 32, 1)
u = 0.5 * sin.(2*pi*85/8192 * (0:8191))
v, = simulateDSM(u, H)

t = 0:100
step(t, v[t.+1])
step(t, u[t.+1])
```
![plot](/example.png)

### Demos

The toolbox also contains a collection of demos as is present in the MATLAB
toolbox. These can be run using `dsdemo#()`. Currently the demos available are:
* `dsdemo1()` - noise transfer function synthesis
* `dsdemo2()` - time-domain simulation, SNR prediction, spectral analysis and SNR calculation
* `dsdemo3()` - coefficient calculation and dynamic range scaling
* `dsdemo5()` - simulation of the element selection logic of a mismatch-shaping DAC
* `dsdemo6()` - hardware-efficient halfband filter design and simulation (there is currently an issue with this, see [#6](/../../issues/6))
* `dsdemo7()` - positively-invariant set computation
