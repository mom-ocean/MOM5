[![CI](https://github.com/mom-ocean/FMS/actions/workflows/CI.yml/badge.svg)](https://github.com/mom-ocean/FMS/actions/workflows/CI.yml)

# Flexible Modeling System (FMS)

The Flexible Modeling System (FMS) is a software framework for supporting the 
efficient development, construction, execution, and scientific interpretation
of atmospheric, oceanic, and climate system models.

More information is available on the [GFDL FMS page](http://www.gfdl.noaa.gov/fms).

# MOM5 Version

This is a [GitHub fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) of 
the [Official FMS repository](https://github.com/NOAA-GFDL/FMS) for the purpose of maintaining
a version of FMS that is compatible with the MOM5 ocean model.

The GFDL `master` branch is tracked in the `https://github.com/mom-ocean/FMS/tree/gfdl/master` branch
on this repository.

The `master` branch in this repository is included as a [subtree](https://www.atlassian.com/git/tutorials/git-subtree)
within the MOM5 repository. The forked occurred from the commit for the 
[Warsaw 201803 Release](https://github.com/mom-ocean/FMS/commit/e8940fe90d68c3dc7c0d6bf1b8f552a577251754). After
this point the codebases started to diverge and it was decided it would be too difficult to udpate 
driver code for all the supported model configurations.

Note that this is not the version of FMS that is used for [MOM6](https://github.com/mom-ocean/MOM6).
# Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an 'as is' basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.

This project code is made available through GitHub but is managed by NOAA-GFDL
at https://gitlab.gfdl.noaa.gov.
