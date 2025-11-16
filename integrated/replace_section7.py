#!/usr/bin/env python3
"""Replace TODO sections in integrated-scfa-microbiome.Rmd with full implementation."""

import re

# Read the implementation sections
with open('integrated-scfa-microbiome-section7-complete.txt', 'r') as f:
    replacement_text = f.read()

# Read the main file
with open('integrated-scfa-microbiome.Rmd', 'r') as f:
    content = f.read()

# Find and replace the TODO sections
# Pattern to match from "# 7. Integrated Analysis:" through the end of "# 8. Summary" TODO
pattern = r'```\{r taxon-scfa-integration.*?```\s*```\{r responder-nonresponder.*?```\s*# 8\. Summary\s*```\{r final-summary.*?```'

# Replace with the new implementation
new_content = re.sub(pattern, replacement_text, content, flags=re.DOTALL)

# Write back
with open('integrated-scfa-microbiome.Rmd', 'w') as f:
    f.write(new_content)

print("Replacement complete")

