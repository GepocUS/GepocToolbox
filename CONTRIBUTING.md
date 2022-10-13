# Welcome to the GepocToolbox contributing guidelines

This document will guide you through the basic steps and the requirements for contributing to this toolbox.

An important and useful way to contribute is to participate in the [discussions](https://github.com/GepocUS/GepocToolbox/discussions) and [issues](https://github.com/GepocUS/GepocToolbox/issues) tabs of the GitHub repository.
For questions or comments on these guidelines, please head to its [discussion thread](https://github.com/GepocUS/GepocToolbox/discussions/1).

## How to contribute

This document is not intended as a guide on how to contribute code to a GitHub repository. You will need to install _git_, basic _git_ knowledge, a GitHub account, a fork of this repository and basic knowledge on how to make _pull requests_. We refer the reader to git's [main page](https://git-scm.com/) for how to install and use _git_ and to GitHub's guides on how to [fork repositories](https://docs.github.com/en/get-started/quickstart/fork-a-repo), make [pull requests](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) and [contribute](https://docs.github.com/en/github/collaborating-with-pull-requests) to other repositories.

## Contribution requirements

All contributions must satisfy the naming conventions, code style and all the other requirements we list. Failure to comply with the requirements we list will result in the rejection of the contribution.

#### General requirements

- All code, comment, documentation, variable names, etc. must be written in English.
- Code should always be accompanied by helpful comments and documentation.
- Coding style, naming conventions, etc. of new files should try to be similar to the ones of other similar files of the toolbox. For instance, all QP solvers should be similar in how they are called and what outputs they return, as well as be coded using a similar style and naming conventions.
- Files should be included in the proper folder or subfolder of the toolbox (see the section "Toolbox structure" down below).
- Empty lines and tabulations should be included to provide a clean, visually pleasing and easy to read/interpret code.
- This toolbox is for academic research and should be simple to use. We want functions and method that are easy to use "out of the box" but that can provide lost of information and customization if so desired. See the section "Arguments" down below for more information in this regard.
- Example files showing how to use, interpret and understand functions and classes should be added if possible.

#### Naming conventions

- Function and variable names should start with a lower-case letter (exceptions such as upper-case letters due to the use of acronyms are accepted). If the name is formed by various words, then they should either be separated by underscores or by using upper-case for the first letter of each word. For example, a function that adds two matrices may be called `add_matrices` or `addMatrices`, but not `AddMatrices`.
- Classes, on the other hand, should start with an upper-case letter if possible.
- Names should be kept simple but provide information to the user when possible. That is, avoid cryptic names.
- Counter variables should be assigned simple names such as: `i`, `j`, `k`, etc.
- Names reserved by standard Matlab functions are not allowed.
- Functions or classes only used within the toolbox should be named starting by `gt_` (standing for _GepocToolbox_). This is not required for functions contained in the `+gt_utils` subfolder (see the section "Toolbox structure").

#### Documentation and comments

- All functions should include an initial documentation containing the following (optional in brackets []):
    + `%%` function_name [- short descriptive sentence]
    + Brief explanation of the function.
    + [explanation of the different ways in which the function can be used]
    + INPUTS: followed by a detailed list of inputs, including their default options, if applicable.
    + OUTPUTS: followed by a detailed list of the outputs.
    + This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
- All classes should include an initial documentation containing the following (optional in brackets []):
    + `%%` class_name [- short descriptive sentence]
    + Brief explanation of the class. Indicate if it extends another class.
    + Constructor [or list of constructors]
    + Properties: followed by a detailed list of properties, indicating default values and any other relevant information (extended properties are not required).
    + Methods: followed by a detailed list of the methods and their inputs/outputs (extended methods are not required).
    + This class is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
- Code should be commented throughout. Distinct sections of the code should be divided with the double comment delimiter `%%`.

#### Arguments

Functions and methods should be built to provide flexibility and ease of use. In particular, we take the following into consideration.

- Input arguments: The [inputParser](https://www.mathworks.com/help/matlab/ref/inputparser.html) class should be used to handle inputs to complex functions. This is not required for all functions, only those with lost of optional inputs with default values or with various ways in which they can be used. If an input parser is used then you should consider:
    + Required inputs, using `addRequired`.
    + Optional inputs, using `addOptional`.
    + Name-value parameters, using `addParameter`.
- Output arguments: The first output argument should be the one that is deemed as the most important one for the end-user. All other outputs should be optional and designed to contain additional useful information. A few noteworthy ones are:
    + `e_flag`: An output with this name will be added, if applicable, to provide an indication of the termination status of the function. This variable will be assigned an integer which will be positive if the function terminated successfully and negative otherwise. The documentation will state the meaning of each value of `e_flag`.
    + `Hist`: A structure containing historics or other specific variables of interest. As an example, an iterative solver may output a variable `Hist` containing the historics of its iterates.
- Names of inputs and outputs follow the same conventions listed above for inner-function variables.

#### Toolbox structure

The toolbox is structured in different folders, each containing functions and classes, typically of a similar theme. New files should be included in the most appropriate folder. The folders of the toolbox are:

- The root directory: This folder should be avoided if possible. It should only include files related to the toolbox itself.
- `/+gt_utils`: This folder contains functions that are not meant for public use. That is, they are only used by other functions and/or classes of the toolbox.
- `/classes`: All classes should be contained here.
- `/functions`: All general-purpose functions should be included here.
- `/solvers`: Optimization solvers should be included here.
- `/benchmarks`: Functions for generating system models and benchmarks go here.

If justified, additional folders and subfolders can be created.

