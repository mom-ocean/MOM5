# Build the version component library

set srcList = ( version )

set lib_name = "lib_version"

# Set up the version string

source ./update_version.csh

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
