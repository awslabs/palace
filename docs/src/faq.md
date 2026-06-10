```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Frequently Asked Questions

```@contents
Pages = ["faq.md"]
Depth = 3
```

This page addresses common questions and issues that users encounter when
working with Palace. If you don't find the answer to your question here, please
check the [issue tracker](https://github.com/awslabs/palace/issues) or consider
opening a new issue.

## Installation

### Can I install Palace on Windows?

We do not officially support or routinely test Palace on Windows and do not
provide precompiled binaries for this platform. Users have reported success
running Palace on [Windows Subsystem for Linux
(WSL)](https://learn.microsoft.com/en-us/windows/wsl/) using an Ubuntu image.
This approach essentially gives you a Linux environment within Windows, which is
compatible with Palace's build system and dependencies.

Support for Palace on WSL is provided on a best-effort basis. If you encounter
problems, please open an issue on our [GitHub
repository](https://github.com/awslabs/palace/issues) and we will do our best to
help you resolve them (but we do not guarantee support).

Before reporting issues, we recommend searching through the [issue
tracker](https://github.com/awslabs/palace/issues) for existing discussions,
known issues, and user-contributed solutions.

##### Common WSL Pitfalls

When working with Palace on WSL, be aware of these two frequent issues:

 1. Missing system packages: The Ubuntu base image in WSL is minimal and lacks
    several packages that Palace's build process expects. You may need to install
    additional system packages even when using Spack, which normally handles
    dependencies automatically. Common missing packages include build tools,
    development headers, and system libraries.

 2. Resource limitations: Windows can constrain the memory, CPU, and disk
    resources allocated to WSL, which may limit the size and complexity of
    simulations you can run. If you encounter out-of-memory errors or unusually
    slow performance, you may need to adjust these limits through the [WSL
    configuration
    settings](https://learn.microsoft.com/en-us/windows/wsl/wsl-config?source=recommendations#main-wsl-settings).

### I installed Palace and downloaded the example files, but I get errors

Palace is under active development with new features constantly being added and
the interface being refined. At this stage, Palace versions are not fully
cross-compatible, which means configuration files and examples from one version
may not work correctly with another version. For this reason, it is crucial to
consult the documentation that corresponds to the exact version of Palace you
are using.

The [official documentation](https://awslabs.github.io/palace/stable/) includes
a version selector at the bottom of the sidebar to help you choose the correct
documentation version. We generally recommend using the most recent version of
Palace unless you have specific reasons to use an older version.

The most common cause of errors when running examples is a mismatch between your
Palace version and the example files you're using. Here's how to ensure
compatibility:

 1. Check your Palace version by running `palace --version` in your terminal
    (if `--version` is supported, if not, check the log of any Palace simulation)
 2. Navigate to the documentation for your specific version using the version selector
 3. Download example files from the GitHub repository at the tag or commit corresponding to your version

!!! note "Following GitHub links from the released documentation"
    
    Currently, the documentation links always point to the GitHub pages for the
    `main` branch. If you follow these links directly, you might get versions of the
    examples that are inconsistent with your Palace installation. To get the correct
    examples:
    
     1. Go to the [Palace GitHub repository](https://github.com/awslabs/palace)
     2. Use the branch/tag selector to switch to the tag corresponding to your version (e.g., `v0.14.0`)
     3. Download the examples from that specific version

### Spack does not have the most recent version of Palace

If you tried installing Palace with Spack and found that you got an outdated
version, there are typically two reasons for this issue.

The most common cause is that you're using an old version of the Spack package
recipes. Spack maintains a repository of package definitions that needs to be
updated periodically to include new software versions.

To update your local Spack repository, run:

```bash
spack repo update -b develop builtin
```

After updating, you can check what versions of Palace are now available with:

```bash
spack versions palace
```

This command will show you all the Palace versions that Spack knows about,
including any newly added ones.

A less common reason is that the newest version of Palace hasn't been added to
Spack yet. This can happen when a new Palace release is very recent and the
Spack maintainers haven't had time to update the package definition.

You can check if this is the case by:

 1. Looking at the [open pull requests](https://github.com/spack/spack-packages/pulls?q=is%3Apr+is%3Aopen+palace) to the `spack-packages` repository
 2. Searching for a PR that updates Palace to the version you want
 3. If such a PR exists, you can track its progress and wait for it to be merged

Once the PR is merged, run the update command shown above to pull the changes to
your local copy of Spack.

##### Installing Development Versions

If you need the absolute latest version of Palace (including features not yet in
a released version), you can install the development version directly from the
main branch:

```bash
spack install palace@develop
```

Keep in mind that development versions may be less stable and are not
recommended for production use.

### I installed Palace with Spack and a new version was released. How do I update Palace?

Spack does not have a built-in mechanism to update packages in place. In fact,
Spack is designed to allow multiple versions of the same package to coexist
simultaneously, which is useful for reproducibility and testing different
versions.

When a new version of Palace is released, you'll need to follow these steps to
get the updated version:

###### Step 1: Update Your Spack Repository

First, update your local copy of the Spack package definitions to include the new Palace version:

```bash
spack repo update -b develop builtin
```

###### Step 2: Verify the New Version is Available

Check that the new version is now available in your Spack installation:

```bash
spack versions palace
```

This should show the new version in the list. If it doesn't appear, the new
version may not have been added to Spack yet (see the previous FAQ entry for
more details).

###### Step 3: Remove the Old Version (Optional)

If you do not have a reason to keep the older version, you can remove it first:

```bash
spack uninstall palace
```

Note that if you have multiple versions installed, Spack will ask you to specify
which one to remove. You can also use `spack find palace` to see all installed
versions first.

###### Step 4: Install the New Version

Finally, install the new version of Palace:

```bash
spack install palace
```

By default, this will install the newest available version. If you want to
install a specific version, you can specify it explicitly:

```bash
spack install palace@0.14.0
```

###### Keeping Multiple Versions

If you prefer to keep both the old and new versions (for example, to maintain
compatibility with existing projects), you can skip step 3 and just install the
new version alongside the old one. You can then use `spack load palace@<version>` to switch between versions as needed.

## I called `spack install palace` but I get an error saying that I need at least one solver

If you followed the instructions to install *Palace* with Spack and found an
error like the following,

```
==> Error: failed to concretize palace for the following reasons:
   1. palace: Need at least one sparse direct solver
   2. palace: At least one eigenvalue solver is required
   3. palace: Need at least one sparse direct solver
    required because conflict is triggered when ~mumps~strumpack~superlu-dist 
     required because palace requested explicitly 
    required because conflict constraint
     required because palace requested explicitly 
   4. palace: At least one eigenvalue solver is required
    required because conflict is triggered when ~arpack~slepc 
     required because palace requested explicitly 
    required because conflict constraint
     required because palace requested explicitly
```

chances are that you do not have a Fortran compiler. This is a common occurrence
on Mac's, which does not ship with a Fortran compiler by default.

Some of the dependencies of *Palace* require a Fortran compiler, and Spack
disables such dependencies when a suitable compiler is not found, leading to the
error message posted above.

To check that this is the case, call

```
spack compiler list
```

This lists the compiler that Spack located. Then, call

```
spack compiler info TOOLCHAIN_NAME
```

where `TOOLCHAIN_NAME` is the name of the compiler listed in the output of the
previous command (e.g., `spack compiler info gcc`). This will provide extra information
on the available compilers, for example:

```
[e]  apple-clang@=17.0.0 build_system=bundle platform=darwin os=tahoe target=aarch64

  prefix: /usr
  compilers:
    cc: /usr/bin/clang
    cxx: /usr/bin/clang++
    fortran: None
```

You have a Fortran compiler if `fortran` is mentioned in the output of `spack compiler info` with output that is not `None`.

If you do not have a Fortran compiler, you can install `gfortran` with your
package manager (e.g., `apt install gfortran` or `brew install gfortran`). Then,
call `spack compiler find` to update the list of compilers that Spack is aware
of.
