# Documentation Standards and Guidelines

## Purpose

This document defines the standards for creating and maintaining documentation in the dronetanalysis project. It ensures consistency, scalability, and ease of navigation across all documentation.

---

## Naming Convention

**Format:** `[CATEGORY]-[NUMBER]-[Descriptive-Name].md`

**Categories:**
- `GUIDE` - User-facing workflows, tutorials, how-to guides
- `OPTIMIZATION` - Performance, memory, and resource optimization
- `FIX` - Critical bug fixes and solutions
- `REFERENCE` - Technical deep-dives, statistical methods, methodology
- `RULES` - Meta-documentation (this file)

**Numbering:**
- Each category has **independent numbering** (01, 02, 03...)
- Numbers are only unique within a category
- Numbers indicate recommended reading order within category
- Use `00-RULES-*` prefix for meta-documentation

**Examples:**
- `GUIDE-01-Complete-Workflow.md` - Main workflow guide
- `GUIDE-02-Network-Metrics.md` - Metrics explanation
- `OPTIMIZATION-01-Memory.md` - Memory optimization
- `FIX-01-Critical-Issues-Summary.md` - Fix summary

**Why This Works:**
- ✅ Add `GUIDE-03-*.md` without renumbering OPTIMIZATION files
- ✅ All GUIDEs group together alphabetically
- ✅ Clear category boundaries in file listings
- ✅ Scalable to 10+ docs per category

---

## When to Create Each Document Type

### GUIDE Documents

**Create when:** Explaining how to use the pipeline or interpret results
**Audience:** End users, researchers, analysts
**Examples:**
- Complete workflow from raw data to results
- How to interpret network metrics
- Troubleshooting common issues

**Template:**
```markdown
# [Title]

## Overview
Brief description of what this guide covers.

## Prerequisites
What the user needs before starting.

## Step-by-Step Instructions
Detailed walkthrough with examples.

## Troubleshooting
Common issues and solutions.

## Related Reading
Links to other relevant docs.
```

---

### OPTIMIZATION Documents

**Create when:** Explaining performance tuning or resource management
**Audience:** Users with large datasets or limited resources
**Examples:**
- Memory optimization strategies
- Disk space reduction techniques
- Cluster resource allocation

**Template:**
```markdown
# [Title] Optimization Guide

## Problem
What resource constraint this addresses.

## Solution
Available optimization strategies.

## Usage
How to enable and configure.

## Performance Impact
Benchmarks and trade-offs.

## Related Reading
Links to other optimization docs.
```

---

### FIX Documents

**Create when:** Documenting critical bug fixes or solutions
**Audience:** Users experiencing errors, developers
**Examples:**
- HDF5 attribute size limit fix
- Memory allocation error fix
- Critical pipeline fixes summary

**Template:**
```markdown
# [Title] Fix

## Problem
Exact error message or behavior.

## Root Cause
Technical explanation of why it happens.

## Solution
Code changes and fixes applied.

## Verification
How to test the fix works.

## Related Reading
Links to related fixes or guides.
```

---

### REFERENCE Documents

**Create when:** Providing deep technical documentation
**Audience:** Advanced users, researchers, developers
**Examples:**
- Statistical methodology
- Pipeline architecture
- Algorithm details

**Template:**
```markdown
# [Title] Reference

## Overview
High-level summary.

## Technical Details
In-depth explanation with code examples.

## Implementation
How it's implemented in the pipeline.

## References
Academic papers, external resources.

## Related Reading
Links to practical guides.
```

---

## Cross-Reference Guidelines

1. **Every document must have a "Related Reading" section**
   - Place at the end of the document
   - Include 3-5 most relevant related docs
   - Provide brief description after each link

2. **Use relative links:**
   ```markdown
   [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md)
   ```

3. **Provide context:**
   ```markdown
   - [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline workflow guide
   ```

4. **Bidirectional links:**
   - If Document A links to Document B, Document B should link back to A
   - Creates a web of interconnected documentation

5. **Update README.md:**
   - Add all new docs to the appropriate category table
   - Update troubleshooting index if applicable

---

## Adding New Documentation

### Step 1: Determine Category

Choose the appropriate category based on the document's purpose:
- **GUIDE** - How to use or do something
- **OPTIMIZATION** - How to improve performance/efficiency
- **FIX** - How to solve a specific problem or error
- **REFERENCE** - Deep technical information

### Step 2: Assign Number

Use the next available number in that category:

**Current documentation:**
```
GUIDE-01-Complete-Workflow.md
GUIDE-02-Network-Metrics.md
OPTIMIZATION-01-Memory.md
OPTIMIZATION-02-Storage.md
OPTIMIZATION-03-Network-Reconstruction.md
FIX-01-Critical-Issues-Summary.md
FIX-02-HDF5-Attributes.md
REFERENCE-01-Statistical-Methods.md
```

**Examples:**
- Adding a new GUIDE → `GUIDE-03-[Name].md`
- Adding a new OPTIMIZATION → `OPTIMIZATION-04-[Name].md`
- Adding a new FIX → `FIX-03-[Name].md`

### Step 3: Use Template

Copy the appropriate template from this document (see "When to Create Each Document Type" section above).

### Step 4: Write Content

Follow the template structure and include:
- Clear, concise explanations
- Code examples where applicable
- Screenshots or diagrams if helpful
- Command-line examples with expected output

### Step 5: Add Cross-References

In the "Related Reading" section, link to:
- The main workflow guide (if not already a guide)
- Related optimization/fix documents
- Technical references if applicable

### Step 6: Update README.md

1. Open [README.md](README.md)
2. Find the appropriate category table
3. Add your document with a brief description
4. Update "Recent Updates" section if it's a major addition

### Step 7: Verify Links

- Check that all internal links work
- Verify cross-references are bidirectional
- Ensure README.md link works

---

## README.md Maintenance

The [README.md](README.md) serves as the central navigation hub. Keep it updated:

### When Adding New Documentation

1. **Update category table:**
   ```markdown
   | Document | Description |
   |----------|-------------|
   | [GUIDE-03-New-Guide.md](GUIDE-03-New-Guide.md) | Brief description |
   ```

2. **Update troubleshooting index** (if applicable):
   ```markdown
   | Error | Document | Section |
   |-------|----------|---------|
   | New error pattern | [FIX-03-New-Fix.md](FIX-03-New-Fix.md) | Problem |
   ```

3. **Update "Recent Updates" section:**
   ```markdown
   ### 2026-02-XX
   - ✅ **New Feature**: Description of what was added
   ```

### When Renaming Files

1. Search and replace all occurrences of the old filename
2. Check all tables and lists
3. Verify all inline links
4. Test navigation from README

### When Archiving Documents

1. Remove from active category tables
2. Add to "Archived Documents" section
3. Note why it was archived (if not obvious)

---

## Archived Documents

Development logs and planning documents go in `archive/`:

**Naming Convention:**
- `DEV-LOG-*` for development logs
- `DEV-PLAN-*` for planning documents

**Purpose:**
- Historical records, not active documentation
- Useful for understanding past decisions
- Reference for similar future work

**Maintenance:**
- Don't update cross-references to archived docs
- Don't include in README.md category tables
- List in "Archived Documents" section of README

---

## Best Practices

### Writing Style

1. **Be concise** - Users want quick answers
2. **Use examples** - Show, don't just tell
3. **Be specific** - "Use `--batch-size 100`" not "use batch mode"
4. **Test your examples** - Ensure code snippets work
5. **Update dates** - Add "Last Updated" at the bottom

### Organization

1. **Start with overview** - Brief summary of what the doc covers
2. **Use clear headers** - Descriptive section titles
3. **Include table of contents** - For docs over 200 lines
4. **End with cross-references** - "Related Reading" section

### Maintenance

1. **Review quarterly** - Check for outdated information
2. **Update when code changes** - Keep examples current
3. **Fix broken links** - Verify links after renaming files
4. **Archive old docs** - Move obsolete docs to archive/

---

## Document Lifecycle

### Creating
1. Choose category and number
2. Use template
3. Write content
4. Add cross-references
5. Update README.md

### Updating
1. Add "Last Updated" date at bottom
2. Update "Recent Updates" in README.md if major
3. Verify examples still work
4. Check that links still valid

### Archiving
1. Move to `archive/` directory
2. Update README.md (remove from active, add to archived)
3. Remove cross-references from active docs
4. Add note explaining why archived

---

## Example: Adding a New GUIDE

**Scenario:** Want to add "GUIDE-03-Cluster-Setup.md" for SGE cluster configuration.

**Steps:**

1. **Determine category:** GUIDE (explains how to do something)

2. **Assign number:**
   ```bash
   ls GUIDE-*.md
   # Shows: GUIDE-01-Complete-Workflow.md, GUIDE-02-Network-Metrics.md
   # Next number: 03
   ```

3. **Create file:** `GUIDE-03-Cluster-Setup.md`

4. **Use template:**
   ```markdown
   # SGE Cluster Setup Guide

   ## Overview
   This guide explains how to configure an SGE cluster for running the dronetanalysis pipeline.

   ## Prerequisites
   - Access to SGE cluster with qsub/qstat commands
   - Python 3.8+ installed on all cluster nodes
   - Shared filesystem accessible from all nodes

   ## Step-by-Step Instructions
   ...

   ## Related Reading
   - [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline workflow
   - [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Memory allocation settings
   ```

5. **Update README.md:**
   ```markdown
   ### 📘 User Guides

   | Document | Description |
   |----------|-------------|
   | [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) | Complete pipeline workflow |
   | [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) | Network topology metrics |
   | [GUIDE-03-Cluster-Setup.md](GUIDE-03-Cluster-Setup.md) | SGE cluster configuration |
   ```

6. **Add cross-references:**
   - In GUIDE-01-Complete-Workflow.md, add link to GUIDE-03
   - In OPTIMIZATION-01-Memory.md, add link to GUIDE-03

7. **Verify:**
   ```bash
   # Check all links work
   grep -r "GUIDE-03" docs/
   ```

**Result:** New guide is integrated, discoverable, and cross-referenced!

---

## Quick Reference

### Current Documentation Structure

```
docs/
├── 00-RULES-Documentation-Standards.md  ← You are here
├── GUIDE-01-Complete-Workflow.md
├── GUIDE-02-Network-Metrics.md
├── OPTIMIZATION-01-Memory.md
├── OPTIMIZATION-02-Storage.md
├── OPTIMIZATION-03-Network-Reconstruction.md
├── FIX-01-Critical-Issues-Summary.md
├── FIX-02-HDF5-Attributes.md
├── REFERENCE-01-Statistical-Methods.md
├── README.md
└── archive/
    ├── DEV-LOG-Bootstrap-Implementation.md
    ├── DEV-LOG-Rewiring-Implementation.md
    └── DEV-PLAN-Focus-Gene-Collection.md
```

### Next Available Numbers

- **GUIDE:** 03
- **OPTIMIZATION:** 04
- **FIX:** 03
- **REFERENCE:** 02

*(Update this section when adding new documents)*

---

## Questions?

If you're unsure about:
- **Which category** - Check "When to Create Each Document Type"
- **How to structure** - Use the appropriate template
- **Cross-references** - Link to related docs in same or different categories
- **README updates** - Follow "README.md Maintenance" section

When in doubt, look at existing documents as examples!

---

**Last Updated:** 2026-02-13
**Status:** ✅ Active documentation standard
**Related Reading:**
- [README.md](README.md) - Central navigation hub
- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Example of GUIDE format
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Example of OPTIMIZATION format
- [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) - Example of FIX format
