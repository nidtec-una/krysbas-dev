# Developing guide for KrySBAS

This documents aims at providing a simple guide to work, develop, and condtribute to KrysBAS. 

- [Create your conda environment](#create-your-conda-environment)
- [Installing required packages](#installing-required-packages)
- [Testing locally with MOxUnit](#testing-locally-with-moxunit)
- [Tagging your commits](#taggin-your-commits)
- [Opening a pull request](#opening-a-pull-request)

## Create your conda environment

To start developing locally, we highly encourage you to work on a dedicated `conda` environment (`venv` + `pip` is also a good alternative you can check). To install `conda`, you can go to [anaconda.org](https://www.anaconda.com/download) and donwload the installer.

Once you downloaded the installer and installed `conda`, you can check whether the installation was succesfull by typing in your terminal 

```
conda list
```

This should show a list of all installed Python packages that come with the base `conda` environment.

Make sure your conda channel is updated by typing

```
conda update -n base -c default conda
```

Now, you can create a new enviroment with Python 3.11:

```
conda create --name krysbas_dev python=3.11
```

To use the environment, you first need to activate it by typing

```
conda activate krysbas_dev
```

If you want to deactivate the environment, use `deactivate` instead.

Although not mandatory, we recommend you to install the following Python packages:

```
conda install git jupyter notebook
```

## Installing required packages

To install the required packages, you first need to clone the KrySBAS repository (do not clone the repo in cloud services such as Dropbox, Google Drive, OneDrive, etc. as this often creates conflicts with Git):

```
git clone https://github.com/nidtec-una/krysbas-dev.git
cd krysbas-dev
```

Now, you can install the required packages by typing

```
conda install --file requirements.txt
```

This will install packages that are used for code styling, documentation, and notebooks.

## Testing locally with MOxUnit

KrySBAS relies heavily on a suite of tests (you can find them in the `tests` directory) to minimize the chances of unwanted bugs. If you want to make a new contribution to the repository (or any change for that matter) you have to make sure that all existing tests pass and provide new tests if you want to add new functionality (reach out to one of the code owners if this is the case).

We employ the [MOxUnit Test Framework](https://github.com/MOxUnit/MOxUnit) to write and run unit tests. Please, follow the [installation instructions](https://github.com/MOxUnit/MOxUnit#installation) to able to run MOxUnit inside MATLAB locally.

## Tagging your commits

To help the code reviewers visualize the changes, we encourage you to tag your commit messages with one of the following tags that best suit the changes that your commit is incorporating:

|Tag|Description|
|----|------------|
|BLD|A commit related to the "build" of the repository, i.e., changes or additions of a new workflow or a required package.|
|BUG|A commit that identifies a bug in the code.|
|BUGFIX|A commit that fixes a bug in the code.|
|ENH|A commit that introduces an enhancement, feature, or new functionality|
|DOC|A commit that adds, fixes, or improves the documenation.|
|MAINT|A commit that introduces a maintainance to the code, i.e., fixing typos, renaming variables, improve logic, etc.|
|STY|A commit related to code style, i.e., changes to please `miss_hit`||
|TST|A commit specifically related to a test.|
|WIP|A commit that is a Work In Progress. Please, be clear on what you are currently working on.|

### Example 1:

Developer X has fixed the docstring of `pd_gmres.m` and writes `DOC: Add description for rel_res output variable`.

### Example 2:

Developer Y has forgotten to make a fix related to the number of lines in a row which caused a `miss_hit` error and writes: `STY: Fix exceeded number of lines`.

## Opening a pull request

Before opening, we kindly ask you to check whether all tests pass and that there are no code style-related issues.

1. From inside MATLAB, go to the `tests` folder and write in the command window
```
run_moxunit_tests
```
This command will run the complete test suite and will let you know whether all tests are succesfull or not.

2. From your terminal (we are assuming you are in the `krysbas_dev` environment), go to the root folder, and type
```
mh_style
```
This command will check whether the M files are compliant with the code style policy. If you see any errors, please go ahead and fix them.

3. Only if you have completed Steps 1 and 2 open a PR. Please, add a succint description on the changes the PR is introducing to the repository, and then a list of specific changes.


