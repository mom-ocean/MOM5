#!/usr/bin/env csh

set root          = $cwd:h                            # The directory you created when you checkout
set code_dir      = $root/src                         # source code directory

# Set up the version string

# The version string is the git hash of the commit used to build the code.
# The version of an executable can be found with the following command:
# strings <executable> | grep 'MOM_COMMIT_HASH='
setenv GIT_CONFIG_NOGLOBAL 'yes'

set old_hash=`grep 'public :: MOM_COMMIT_HASH' ../src/version/version.F90 | cut -d '"' -f 2 | cut -d '=' -f 2`
set new_hash=`git rev-parse HEAD`

if ( $old_hash != $new_hash ) then
    echo "Current version hash:  $old_hash"
    echo "Updating version hash: $new_hash"
    sed -e "s/{MOM_COMMIT_HASH}/$new_hash/g" $code_dir/version/version.F90.template > $code_dir/version/version.F90
endif

