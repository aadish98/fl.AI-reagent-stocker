#!/usr/bin/env node
/** Build priority-run overview deck via python-pptx (PowerPoint-native OOXML). */
const path = require("path");
const { execFileSync } = require("child_process");

const RUN = "priority_run_06-30-2026_v3";
const out = path.join(__dirname, "..", "data", "gene_sets", RUN, "Priority_Gene_Sets_Overview.pptx");

execFileSync("python3", [path.join(__dirname, "build_priority_run_v3_pptx.py")], { stdio: "inherit" });
console.log("WROTE " + out);
