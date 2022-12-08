# Mesh Converter

## Running

```bash
./converter path [type] [options]
```

## Gettting help

```bash
./converter --help
```

## Building

If the code is changed run `make build` to re-build the docker image

## Exemples

### Malha Petrobras

```bash
$ ./converter malhas/malha_petrobras/MESH.txt 2 --wells-path malhas/malha_petrobras/INPUT.lua --inactive-cells-path malhas/malha_petrobras/ACTNUM.txt --tetra-optimization --refinement-threshold 0.150
```

### Malha Namorado

```bash
$ ./converter malhas/malha_namorados/MESH3.txt 2 --refinament-threshold 0.01 --deactivate-hexa-collapsing --tetra-optimization
```
