```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Working with Spack

[Spack](https://spack.io/) is a multi-platform package manager designed to
handle the complexity of building scientific software with multiple versions,
configurations, and dependencies. While Spack is commonly used to install
*Palace* as an end user would, it's also a powerful tool for active development.
One of the main reasons you might consider doing so is that Spack makes
co-development of multiple packages straightforward, e.g., allowing you to point
to local copies of both *Palace* and MFEM simultaneously. This guide walks you
through both using Spack to develop *Palace* itself and working on *Palace*'s
Spack package recipe.

Before diving in, this guide assumes you already have Spack installed and
understand the basics of how it works. If you're completely new to Spack, you
should start with their official documentation to get familiar with concepts
like specs, variants, and environments. Check the [official
tutorials](https://spack-tutorial.readthedocs.io/en/latest/index.html) out to
get started.

## Developing Palace with Spack Environments

The recommended approach for developing *Palace* with Spack is to use Spack
environments. An environment is essentially a self-contained workspace defined
by a `spack.yaml` file that specifies exactly what you want to build and how you
want to build it. This keeps your development setup isolated and reproducible.

To get started, create a dedicated directory for your environment. You might
call it something like `palace_spack` or `palace_dev`. Inside this directory,
create a `spack.yaml` file with the following content:

```yaml
spack:
  specs:
  - palace@develop
  develop:
    palace:
      spec: palace@=develop
      path: ~/repos/palace
```

The `path` field should point to wherever you've cloned the *Palace* repository
on your system. The `specs` section defines what you want to install, and you
can customize this with variants. For instance, if you need the test suite
built, you would change the spec to `palace@develop+tests`. Spack's variant
system is quite flexible, so you can enable or disable features as needed for
your development work.

Before running the installation, it's sometimes worth taking advantage of
packages that are already installed on your system. Spack calls these "external
packages," and finding them can significantly speed up your build times. You can
scan for them by running:

```sh
spack -e palace_spack external find --all
```

This command will modify your `spack.yaml` file to include references to system
packages like compilers, MPI implementations, and common libraries. Spack will
then use these instead of building them from scratch, which can save hours on
the first build.

!!! note "Risks of using externals"
    
    While using external packages can save hours of compilation time, Spack
    does not always handle them correctly, and incompatibilities are possible.
    If you find issues, you can try removing the external packages and letting
    Spack compile everythin.

Now you're ready to install *Palace* by running:

```sh
spack -e palace_spack install 
```

Be prepared for this first installation to take some time, potentially a few
hours depending on your system and what needs to be built. Spack is building not
just *Palace* but all of its dependencies from source. Spack caches everything
it builds, so subsequent installations will be much faster. Even if you change
variants or make other modifications to your environment, Spack is smart enough
to only rebuild what's actually affected by your changes.

Once the installation completes, you can activate the environment and load
*Palace* into your shell:

```sh
spack env activate palace_spack
spack load palace
```

After this, the `palace` command will be available in your terminal. Whenever
you make changes to the *Palace* source code, you can rebuild by simply running
`spack install` again from within the environment. Spack will detect what's
changed and recompile only what's necessary.

The beauty of this development mode is that it extends to other packages as
well. If you need to modify MFEM alongside *Palace*, you can add another entry
to the `develop` section of your `spack.yaml` pointing to your local MFEM
checkout. Spack will then manage both packages in development mode, rebuilding
each as needed when you make changes.

## Developing Palace's Package Recipe

Beyond developing *Palace* itself, you might need to work on *Palace*'s Spack
package recipe (the `package.py` file that tells Spack how to build *Palace*).
This is particularly relevant if you're adding new dependencies, changing build
options, or preparing a new release.

The first step is to add *Palace*'s Spack repository to your Spack
configuration. Assuming your *Palace* checkout is in `~/repos/palace`, you would
run:

```sh
spack repo add ~/repos/palace/spack_repo/local
```

You should see a confirmation message like `==> Added repo with namespace 'local'.` The namespace `local` is significant because it allows you to
distinguish between your development version of the package and the official one
in Spack's built-in repository. When you refer to `local.palace`, you're
explicitly using the recipe from your personal repository, while
`builtin.palace` refers to the official recipe that ships with Spack.

To verify that everything is working, try making a small change to the package
file. For example, you could edit `spack_repo/local/packages/palace/package.py`
and change the docstring to all capital letters. Then run:

```sh
spack info local.palace | head -n 4
```

You should see your modified docstring. Compare this with the output of `spack info builtin.palace | head -n 4` to see the difference between your local
version and the official one.

A shortcut to end the `package.py` is calling `spack edit local.palace`, which
will open the file with your `$EDITOR`.

One important detail is that the `local` repository takes priority over
`builtin`. This means that if you just run `spack info palace` without a prefix,
Spack will use your local version. While this is convenient, we'll continue
using the explicit `local.` prefix in examples for clarity.

#### Working Around a Spack Bug with Patches

There's currently a known issue in Spack version <=1.2.0 that affects packages
using patches when they come from non-builtin repositories (see the [bug
report](https://github.com/spack/spack/issues/51505)). Since *Palace* uses
patches, you'll run into this problem. The workaround is straightforward but a
bit inelegant: you need to add a copy of the MFEM package to your local
repository as well.

Rather than maintaining a duplicate copy of the MFEM recipe, you can create a
symbolic link from Spack's built-in MFEM package to your local repository:

```sh
ln -s "$(spack location --repo builtin)"/packages/mfem ~/repos/palace/spack_repo/local/packages/
```

This command uses Spack's `location` subcommand to find where the built-in
repository lives, then creates a symlink to the MFEM package directory. To
verify that the symlink was created correctly, run:

```sh
spack info local.mfem | head -n 4
```

If you see MFEM's package information, you're all set. This workaround ensures
that when Spack processes patches for *Palace*, it can find MFEM in the same
repository namespace and avoid the bug.

## Iterating on Package Changes

With your local repository configured, you can now iterate on *Palace*'s
`package.py` file. After making changes, you can test them by installing
*Palace* again. You can also use `spack spec local.palace@develop` to see how
Spack would resolve dependencies without actually building anything, which is
useful for catching errors early.

### Working with source code patches

*Palace* requires some patches in some of the dependencies (primarily MFEM).
Spack handles this, but there are a few common pitfalls.

First, patches are evaluated at concretization. This is because patches are
embedded in the hash, which is computed at concretization. If you are iterating
on patches, make sure to re-concretize your environment.

Second, *Palace* depends on a specific commit of MFEM (defined in the
`package.py` or in the `ExternalGitTags.cmake`). `spack develop mfem` does not
respect this commit restriction. In most cases, this means that you have to
manually git checkout the specific hash for mfem in the mfem folder so that
patches apply cleanly.

### Working with patches for `package.py` files

When working with Spack, it sometimes happens that we find problems in upstream
`package.py` files. When this happens, we open pull requests to
[`spack-packages`](https://github.com/spack/spack-packages) to fix the issue.
However, having those pull requests merged can take weeks. To ensure that our
continuous integration works, we manually apply those PR as part setting Spack
up in our runners (the `setup-runner` workflow).

If you are working with Spack and *Palace* locally, you might need to apply the
same patches. To apply a patch in `spack_repo/patches`:

```bash
cd $(spack location --repo builtin)
git apply /path/to/palace/spack_repo/patches/prNUM.diff
```

Note that if you do so, your `builtin` repository will no longer
be clean and you will not be able to update `spack-packages` as easily.

Instead, if you some patches applied and need to upgrade the recipes, the
easiest way is probably to remove all the patches and re-apply the ones you need:

```bash
cd $(spack location --repo builtin)
git fetch
git reset --hard origin/develop
# Apply relevant patches
```

Patches are often more less common build options, so you might not need all of
them.
