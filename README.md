## Examples commands to run multiz

```bash
multiz human.chimp.galago.maf human.mouse.rat.maf 1 
```

## program descriptions

```text
multiz.v11.2:  -- aligning two files of alignment blocks where top rows are always the reference, reference in both files cannot have duplicats
args: [R=?] [M=?] file1 file2 v? [out1 out2] [nohead] [all]
	R(30) radius in dynamic programming.
	M(1) minimum output width.
	out1 out2(null) null: stdout; out1 out2: file names for collecting unused input.
	nohead(null) null: output maf header; nohead: not to output maf header.
	all(null) null: not to output single-row blocks; all: output all blocks.

maf-file1 and maf-file2 are two maf files to be aligned, each 
topped by a same reference sequence. The alignment of reference 
sequence with other components might be just for purpose of 
determing approximate alignment between two files, thus the 
alignment might be fixed or not, this is specified by v value, 
which can be only 0 or 1.

0 - neither alignment of reference in each file is fixed.
1 - the alignment of reference in the first file is fixed.

[R=?] species radius values in dynamic programming, by default 30
[M=?] species minimum output width, by default 1, which means 
output all blocks.

[out1] collects unused blocks from maf-file1
[out2] collects unused blocks from maf-file2
[nohead] specifies not to have maf header for output
```
