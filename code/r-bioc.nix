#Copyright 2020 RICHARD TJÖRNHAMMAR
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.


# https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/r-modules/bioc-packages.nix

let
  pkgs = import <nixpkgs> {};
  stdenv = pkgs.stdenv;
in with pkgs; {
  myProject = stdenv.mkDerivation {
    name = "Speedracer";
    version = "1";
    buildInputs =  [
      R
      rPackages .knitr
      rPackages .rmarkdown
      rPackages .oligo
      rPackages .pd_hta_2_0
      rPackages .DOSE
      rPackages .Cairo
      rPackages .org_Hs_eg_db
      rPackages .clusterProfiler
    ];
  };
}

