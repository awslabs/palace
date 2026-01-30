# PR Stack Code Reviews

Reviews for the config refactoring PR stack, PRs 2-6.

---

## PR 2: `hughcars/remove-redundant-config-checks`

**Stats:** 7 files, +151/-1042 lines (net -891 lines)

### Summary
Removes redundant runtime validation checks from `configfile.cpp` that are now enforced by JSON schema validation. Also adds if/then schema rules for Problem.Type → Solver requirements.

### Strengths
1. **Significant code reduction** - Removes ~900 lines of redundant validation code
2. **Schema improvements** - Adds proper numeric constraints (`exclusiveMinimum: 0` for Index fields, `minimum/maximum` for Order)
3. **if/then rules** - Properly enforces Problem.Type → Solver section requirements at schema level
4. **Consistent API change** - Changes `SetUp(json&)` to `SetUp(const json&)` throughout, removing mutation

### Issues/Concerns
1. **Incomplete removal of JSON_DEBUG** - The `constexpr bool JSON_DEBUG = false;` is removed but some debug blocks may remain in later code
2. **Error message quality** - Schema validation errors may be less user-friendly than the original MFEM_VERIFY messages (e.g., "instance is below or equals minimum of 0" vs "Index should be an integer > 0")
3. **Missing schema constraint** - `Absorbing.Order` is constrained to `[1,2]` but the original code had a more descriptive error message

### Suggestions
1. Consider adding `"errorMessage"` annotations to schema for better user-facing errors (if the validator supports it)
2. Verify all `JSON_DEBUG` blocks are removed (search for any remaining)
3. Add schema test for the if/then Problem.Type → Solver rules

### Verdict: **Approve**
The PR achieves its goal of removing redundant checks. The schema constraints are correct and the code is cleaner.

---

## PR 3: `hughcars/refactor-config-structs`

**Stats:** 3 files, +246/-373 lines (net -127 lines)

### Summary
Removes `DataVector` and `DataMap` template abstractions, replacing them with direct use of `std::vector` and `std::map`. Also separates `Terminal` from `LumpedPort` as a first-class boundary type.

### Strengths
1. **Removes unnecessary abstraction** - `DataVector` and `DataMap` were thin wrappers that added no value
2. **Direct STL usage** - Using `std::vector<MaterialData>` and `std::map<int, ProbeData>` directly is clearer
3. **Modern C++** - Uses structured bindings (`auto [iter, inserted] = ...`) and `try_emplace`
4. **Terminal separation** - Makes the config structure match the actual domain model better

### Issues/Concerns
1. **Inconsistent constructor patterns** - Some structs get constructors (`MaterialData(const json&)`), others don't yet
2. **SetUp methods remain** - Many `SetUp` methods still exist alongside new constructors, creating inconsistency

### Suggestions
1. This PR is a stepping stone - the inconsistency is addressed in PRs 4-5
2. Consider adding a brief comment explaining why Terminal is now separate from LumpedPort

### Verdict: **Approve**
Clean removal of unnecessary abstractions. The inconsistency with SetUp methods is intentionally addressed in subsequent PRs.

---

## PR 4: `hughcars/config-constructor-conversions-boundaries`

**Stats:** 3 files, +196/-315 lines (net -119 lines)

### Summary
Converts boundary-related `SetUp` methods to constructors: `LumpedPortData`, `WavePortData`, `SurfaceCurrentData`, and other boundary data types.

### Strengths
1. **Consistent pattern** - All boundary data types now use constructors
2. **Cleaner initialization** - Objects are fully initialized at construction time
3. **Test updates** - Tests properly updated to use new constructor pattern

### Issues/Concerns
1. **Test verbosity** - Tests now require `config.find("...")->find("Boundaries")` pattern which is more verbose
2. **No default constructors removed** - Default constructors still exist, allowing partially-initialized objects

### Suggestions
1. Consider whether default constructors should be deleted for types that require JSON input
2. The test pattern could be simplified with a helper function

### Verdict: **Approve**
Consistent conversion to constructor pattern. The verbosity in tests is acceptable.

---

## PR 5: `hughcars/config-constructor-conversions-solvers`

**Stats:** 6 files, +substantial changes

### Summary
Converts remaining `SetUp` methods to constructors: solver data types (`DrivenSolverData`, `EigenSolverData`, `LinearSolverData`, etc.), `DomainData`, `BoundaryData`, `ModelData`, `ProblemData`.

### Strengths
1. **Completes the refactor** - All config data types now use constructors
2. **IoData simplified** - Top-level parsing in `IoData` is now cleaner with direct constructor calls
3. **Consistent API** - All config structs follow the same pattern

### Issues/Concerns
1. **Cross-section validation** - The `Eigenmode.Target` requirement when `Driven` section exists is now checked in `SolverData` constructor, which has access to both sections. This is correct but the logic is subtle.
2. **Large PR** - This PR touches many files and could be harder to review in isolation

### Suggestions
1. Add a comment explaining the `Eigenmode.Target` / `Driven` interaction
2. Consider if any cross-field validations should move to a separate `Validate()` function (addressed in PR 6)

### Verdict: **Approve**
Completes the constructor conversion consistently.

---

## PR 6: `hughcars/increased-schema-tests`

**Stats:** 7 files, +459/-150 lines

### Summary
Adds comprehensive schema tests, moves more validation to schema (mutual exclusion, required fields), and adds a `Validate()` function for cross-field checks that can't be expressed in JSON schema.

### Strengths
1. **Extensive test coverage** - Adds tests for array types, required fields, mutual exclusion, error message format
2. **Schema improvements** - PEC/Ground and PMC/ZeroCharge mutual exclusion now in schema using `allOf/not`
3. **oneOf constraints** - LumpedPort and SurfaceCurrent now require either Attributes or Elements
4. **Validate() function** - Clean separation of cross-field validation (duplicate indices, excitation matching)
5. **Bug fix** - Changes `target_upper > target` to `target_upper >= target`

### Issues/Concerns
1. **Removed test** - "Schema Validation - Embedded Schema Matches Source" test was removed. Was this intentional?
2. **Duplicate code** - The excitation validation logic appears in both `BoundaryData` constructor (for normalization) and `Validate()` (for error checking). Could be consolidated.
3. **Validate() not called** - The new `Validate()` function is defined but I don't see where it's called in the main code path

### Suggestions
1. **Critical:** `config::Validate(boundaries)` is defined but NEVER CALLED. Must add call in `IoData` constructor after `boundaries = config::BoundaryData(*boundaries_it);`
2. Restore or explain removal of embedded schema test
3. Consider consolidating excitation validation logic
4. Add test for duplicate index detection across LumpedPort/WavePort/SurfaceCurrent/Terminal

### Additional Schema Opportunities
1. Could add `minItems: 1` to Attributes arrays to catch empty arrays
2. Could add pattern validation for Direction strings (already done via enum, good)
3. Could add `uniqueItems: true` to Attributes arrays

### Verdict: **Request Changes**
The `Validate()` function appears to be defined but not called. This needs to be verified/fixed before merging.

---

## Summary Table

| PR | Branch | Verdict | Key Action Items |
|----|--------|---------|------------------|
| 2 | remove-redundant-config-checks | ✅ Approve | None |
| 3 | refactor-config-structs | ✅ Approve | None |
| 4 | config-constructor-conversions-boundaries | ✅ Approve | None |
| 5 | config-constructor-conversions-solvers | ✅ Approve | Add comment for Target/Driven logic |
| 6 | increased-schema-tests | ⚠️ Request Changes | Verify Validate() is called |

