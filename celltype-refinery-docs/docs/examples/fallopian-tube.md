---
sidebar_position: 1
---

# Fallopian Tube Example

Example annotation workflow for fallopian tube tissue.

## Data

- 3.4M cells
- 45 markers
- 5 donors
- 4 regions (fimbriae, ampulla, isthmus, uterine junction)

## Workflow

1. Preprocessing with batch correction
2. Leiden clustering at resolution 0.6
3. Annotation with FT marker map
4. Iterative refinement (3 rounds)
5. Composition and spatial analysis

## Results

- 91% assignment rate
- 15 final cell types
- Clear regional patterns
