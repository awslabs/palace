```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Configuration File

*Palace* simulations are configured using a single [JSON](https://en.wikipedia.org/wiki/JSON) file
passed as the command-line argument. If you are new to *Palace*, start with the [quick-start
guide](../quick.md) to see a full worked example.

The configuration file has five required top-level sections:

```json
{
    "Problem":    { ... },
    "Model":      { ... },
    "Domains":    { ... },
    "Boundaries": { ... },
    "Solver":     { ... }
}
```

The complete reference for every configuration option is in the [Configuration File
Reference](reference.md). Each entry shows the JSON key, type, default value, and description.
Fields marked **required** have no default and must be supplied.

## File format

The configuration file uses JSON with two custom extensions that make it easier to write by hand.

  - *Comments*: C (`/* */`) and C++ (`// ...`) style comments can be included.
  - *Integer range expansion*: Integer arrays for mesh attributes are comma-separated lists, but can
    include inclusive ranges `a-b`. For example, `"Attributes": [1, 3-5, 8]` is equivalent to
    `"Attributes": [1, 3, 4, 5, 8]`.

These custom extensions are preprocessed into standard JSON before parsing. Apart from these, the
file must be valid JSON. Duplicate keys within any object are detected and reported as an error.

## JSON Validation

*Palace* validates the configuration file against a [JSON
Schema](https://json-schema.org/draft-07/schema) before the solver starts. The same schema can be
used independently to check your configuration before running *Palace* by running:

```bash
./scripts/validate_config config.json
```

Note that there are constraints on the configuration that cannot be expressed in a JSON Schema.
*Palace* may still fail to run with your configuration even if JSON Schema validation passes.
Consult the [Configuration File Reference](reference.md) and the [User Guide](../guide/guide.md) for
the requirements of individual settings.
