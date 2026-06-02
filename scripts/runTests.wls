#!/usr/bin/env wolframscript
(* Runs every Tests/*.wlt suite and exits 1 if any test fails.
   wolframscript exits 0 even on test failure, so the explicit Exit is what
   makes CI go red. Run from anywhere: wolframscript -file scripts/runTests.wls *)

root = DirectoryName[$InputFileName, 2];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];

files = FileNames["*.wlt", FileNameJoin[{root, "Tests"}]];
reports = AssociationMap[TestReport[#] &, files];
KeyValueMap[
  Print[FileNameTake[#1], ": ",
     If[#2["AllTestsSucceeded"], "PASS",
        "FAIL (" <> ToString[#2["TestsFailedCount"]] <> ")"]] &,
  reports];

allOK = AllTrue[Values[reports], #["AllTestsSucceeded"] &];
Print[If[allOK, "ALL PASS", "FAILURES"]];
Exit[If[allOK, 0, 1]]
