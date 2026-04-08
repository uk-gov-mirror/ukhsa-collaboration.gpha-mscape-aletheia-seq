# mscape-template

<<<<<<< HEAD
This repository is a template for creating new repositories containing  
code that will run on mSCAPE. It serves as a guide for code layout  
and files will need amending to fit the repo purpose.  

As a minimum, mSCAPE repositories should include the following:
- Code to be in src/ layout 
- tests/ folder at same level as src/ 
- .github/ folder containing workflows and a pull request template
- .gitignore file
- pyproject.toml file
- .pre-commit-hooks.yaml

This repo follows the above structure and contains examples of the files  
referenced above.

A ssuggested layout for repo READMEs is included below. The guidance  
documentation contains further information on required repository  
structure, development cycles, and making pull requests. Please read  
this guidance document before using the template.

# mscape-template
=======
This repository is a template for creating new repositories containing code that 
will run on mSCAPE. It serves as a guide for code layout and files will need 
amending to fit the repo purpose.

There is guidance documentation including SOPs that contain further information on required repository
structure, development cycles, and making pull requests. Please read
this guidance documentation before using the template.

Below the dashed line is a README template to be edited.

---------------------
# Project or Repo Name
>>>>>>> template/main

Brief description of project here

## Installation

<<<<<<< HEAD
Add installation instructions here. Ideally include commands to make  
the process as easy as possible for users.  

Clone repo and create environment:  
`git clone git@github.com:ukhsa-collaboration/mscape-template.git`  

`conda env create -n mscape_template `  

`conda activate mscape_template`  

Installation for users:  
`cd mscape-template`  
`pip install .`

Installation for developers (installs code in editable mode):  
`cd mscape-template`  
=======
Add installation instructions here. Ideally include commands to make
the process as easy as possible for users.

Clone repo and create environment:

`git clone git@github.com:ukhsa-collaboration/project-name.git`

`conda env create -n project-name `

`conda activate project-name`

Installation for users: 

`cd project-name`

`pip install .`

Installation for developers (installs code in editable mode): 

`cd project-name`

>>>>>>> template/main
`pip install --editable '.[dev]'`

## Usage

<<<<<<< HEAD
Include command line arguments (e.g. the output displayed when using -h)  
for reference. Example commands can also be helpful.

## Other sections

Add other sections as appropriate for your repo. This may include  
instructions on updating the repo, instructions on adding new  
references, troubleshooting etc. 

=======
Include command line arguments (e.g. the output displayed when using -h)
for reference. Example commands can also be helpful.

```
project-name --input <path> --output <path>
```

## Inputs (optional - useful if you have a CLI)

You may wish to use a table to list out the args:

| Argument | Required | Description |
| -------- | ------- | ------- |
| --input, -i | Yes | Input for command line use  |
| --output, -o | Yes | Output for command line use  |

## Outputs (optional - useful if your tool creates lots of output files)

Explain the output files in detail here.

## Other sections

Add other sections as appropriate for your repo. This may include
instructions on updating the repo, instructions on adding new
references, troubleshooting etc.
>>>>>>> template/main
