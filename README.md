[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/D1Eq9Is4)
[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-2e0aaae1b6195c2367325f4f02e2d04e9abb55f0b24a779b69b11b9e10269abc.svg)](https://classroom.github.com/online_ide?assignment_repo_id=21293818&assignment_repo_type=AssignmentRepo)
# pro3. perform global|local alignment
* your name
* student ID

## Description

* Write a Python script to perform a global or local alignment.
* Creating your own script, i.e. hw3.py.
* In this program, library Biostrings is only used to parse input fasta file.
* Packages you can use: numpy, pandas, Bio
* You should write a program with a function named alignment, ie.
```
def alignment(input_path, score_path, output_path, aln, gap):
    .
    .
    .
    .
```
* hw3_ref.py: You can start from this reference code and try to write your own comment in English.
* If there is more than one local alignment with the same highest score, you should output local alignments with the maximum length. 
* If there is more than one local alignment with the same highest score and length, you should output those local alignments in string sequential order according to protein1 and then protein2, i.e., 
  ```
  >protein1
  local alignment1
  >protein2
  local alignment1
  >protein1
  local alignment2
  >protein2
  local alignment2
  ```
## Parameters

* input: .fasta file (ie. examples/test_global.fasta, examples/test_local.fasta)
* score: score file (ie. examples/pam100.txt, examples/pam250.txt)
* aln: global|local
* gap: gap score
* output: .fasta file. You should output your alignment in FASTA format (ie. examples/result_global.fasta, examples/result_local.fasta).
  
## Command

Executing your code with the following command.


```Python
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
```

## Evaluation

10 testing data(5 public, 5 private)

The correct answer gets 10 points for each testing data.


### Penalty

* High code similarity to others: YOUR SCORE = 0

## References
Please provide the code along with its reference. For example, you can cite it as: ```# ChatGPT, respond to “your prompt,” on February 16, 2023```. Below is an example of a reference format summarizing the use of ChatGPT for R programming
提供Gemini 2.5Pro : 
"
write a program with a function named alignment, ie.
def alignment(input_path, score_path, output_path, aln, gap):
If there is more than one local alignment with the same highest score, you should output local alignments with the maximum length.
"並提供input, output 範本
Gemini 提供的python, 針對Local Alignments 部分要調整Smith–Waterman結束的條件，
Gemini 提供在Trace back 遇到矩陣分數0時停止Trace back，並輸出序列，
但在gradescope第三題可看出比對結果為0也屬於結果序列一部分，
所以需調整終止條件，讓比對為0的部分也輸出到結果序列中。





