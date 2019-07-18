# Build the version component library

set srcList = ( version )

set lib_name = "lib_version"

# Set up the version string

# The version string is the git hash of the commit used to build the code.
# The version of an executable can be found with the following command:
# strings <executable> | grep 'MOM_COMMIT_HASH='
setenv GIT_CONFIG_NOGLOBAL 'yes'

set old_hash=`grep 'public :: MOM_COMMIT_HASH' ../src/version/version.F90 | cut -d '"' -f 2 | cut -d '=' -f 2`
set new_hash=`git rev-parse HEAD`

if ( $old_hash != $new_hash ) then
    sed -e "s/{MOM_COMMIT_HASH}/$new_hash/g" $code_dir/version/version.F90.template > $code_dir/version/version.F90
endif

set curdir=$PWD
mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" $srcList
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif

cd $curdir
