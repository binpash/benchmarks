#!/bin/bash

# This function takes no arguments
# It tries to determine the name of this file in a programatic way.
_get_sourced_filename() {
    if [ -n "${BASH_SOURCE[0]}" ]; then
        basename "${BASH_SOURCE[0]}"
    elif [ -n "${(%):-%x}" ]; then
        # in zsh use prompt-style expansion to introspect the same information
        # see http://stackoverflow.com/questions/9901210/bash-source0-equivalent-in-zsh
        basename "${(%):-%x}"
    else
        echo "UNKNOWN FILE"
    fi
}

# The arguments to this are:
# 1. activation nature {activate|deactivate}
# 2. toolchain nature {build|host|ccc}
# 3. machine (should match -dumpmachine)
# 4. prefix (including any final -)
# 5+ program (or environment var comma value)
# The format for 5+ is name{,,value}. If value is specified
#  then name taken to be an environment variable, otherwise
#  it is taken to be a program. In this case, which is used
#  to find the full filename during activation. The original
#  value is stored in environment variable CONDA_BACKUP_NAME
#  For deactivation, the distinction is irrelevant as in all
#  cases NAME simply gets reset to CONDA_BACKUP_NAME.  It is
#  a fatal error if a program is identified but not present.
_tc_activation() {
  local act_nature=$1; shift
  local tc_nature=$1; shift
  local tc_machine=$1; shift
  local tc_prefix=$1; shift
  local thing
  local newval
  local from
  local to
  local pass

  if [ "${act_nature}" = "activate" ]; then
    from=""
    to="CONDA_BACKUP_"
  else
    from="CONDA_BACKUP_"
    to=""
  fi

  for pass in check apply; do
    for thing in $tc_nature,$tc_machine "$@"; do
      case "${thing}" in
        *,*)
          newval=$(echo "${thing}" | sed "s,^[^\,]*\,\(.*\),\1,")
          thing=$(echo "${thing}" | sed "s,^\([^\,]*\)\,.*,\1,")
          ;;
        *)
          newval="${CONDA_PREFIX}/bin/${tc_prefix}${thing}"
          if [ ! -x "${newval}" -a "${pass}" = "check" ]; then
            echo "ERROR: This cross-compiler package contains no program ${newval}"
            return 1
          fi
          ;;
      esac
      if [ "${pass}" = "apply" ]; then
        thing=$(echo ${thing} | tr 'a-z+-' 'A-ZX_')
        eval oldval="\$${from}$thing"
        if [ -n "${oldval}" ]; then
          eval export "${to}'${thing}'=\"${oldval}\""
        else
          eval unset '${to}${thing}'
        fi
        if [ -n "${newval}" ]; then
          eval export "'${from}${thing}=${newval}'"
        else
          eval unset '${from}${thing}'
        fi
      fi
    done
  done
  return 0
}

if [ "${CONDA_BUILD:-0}" = "1" ]; then
  CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem ${PREFIX}/include -fdebug-prefix-map=${SRC_DIR}=/usr/local/src/conda/${PKG_NAME}-${PKG_VERSION} -fdebug-prefix-map=${PREFIX}=/usr/local/src/conda-prefix"
  DEBUG_CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem ${PREFIX}/include -fdebug-prefix-map=${SRC_DIR}=/usr/local/src/conda/${PKG_NAME}-${PKG_VERSION} -fdebug-prefix-map=${PREFIX}=/usr/local/src/conda-prefix"
  LDFLAGS_USED="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,${PREFIX}/lib -L${PREFIX}/lib"
  CPPFLAGS_USED="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem ${PREFIX}/include"
  DEBUG_CPPFLAGS_USED="-D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem ${PREFIX}/include"
  CMAKE_PREFIX_PATH_USED="${PREFIX}:${CONDA_PREFIX}/x86_64-conda_cos6-linux-gnu/sysroot/usr"
else
  CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem ${CONDA_PREFIX}/include"
  DEBUG_CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem ${CONDA_PREFIX}/include"
  CPPFLAGS_USED="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem ${CONDA_PREFIX}/include"
  DEBUG_CPPFLAGS_USED="-D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem ${CONDA_PREFIX}/include"
  LDFLAGS_USED="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,${CONDA_PREFIX}/lib -Wl,-rpath-link,${CONDA_PREFIX}/lib -L${CONDA_PREFIX}/lib"
  CMAKE_PREFIX_PATH_USED="${CONDA_PREFIX}:${CONDA_PREFIX}/x86_64-conda_cos6-linux-gnu/sysroot/usr"
fi

if [ "${CONDA_BUILD:-0}" = "1" ]; then
  if [ -f /tmp/old-env-$$.txt ]; then
    rm -f /tmp/old-env-$$.txt || true
  fi
  env > /tmp/old-env-$$.txt
fi

_tc_activation \
  deactivate host x86_64-conda_cos6-linux-gnu x86_64-conda_cos6-linux-gnu- \
  cc cpp gcc gcc-ar gcc-nm gcc-ranlib \
  "CPPFLAGS,${CPPFLAGS:-${CPPFLAGS_USED}}" \
  "CFLAGS,${CFLAGS:-${CFLAGS_USED}}" \
  "LDFLAGS,${LDFLAGS:-${LDFLAGS_USED}}" \
  "DEBUG_CPPFLAGS,${DEBUG_CPPFLAGS:-${DEBUG_CPPFLAGS_USED}}" \
  "DEBUG_CFLAGS,${DEBUG_CFLAGS:-${DEBUG_CFLAGS_USED}}" \
  "CMAKE_PREFIX_PATH,${CMAKE_PREFIX_PATH_USED}" \
  "_CONDA_PYTHON_SYSCONFIGDATA_NAME,${_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED}" \
  "CONDA_BUILD_SYSROOT,${CONDA_PREFIX}/x86_64-conda_cos6-linux-gnu/sysroot"

if [ $? -ne 0 ]; then
  echo "ERROR: $(_get_sourced_filename) failed, see above for details"
else
  if [ "${CONDA_BUILD:-0}" = "1" ]; then
    if [ -f /tmp/new-env-$$.txt ]; then
      rm -f /tmp/new-env-$$.txt || true
    fi
    env > /tmp/new-env-$$.txt

    echo "INFO: $(_get_sourced_filename) made the following environmental changes:"
    diff -U 0 -rN /tmp/old-env-$$.txt /tmp/new-env-$$.txt | tail -n +4 | grep "^-.*\|^+.*" | grep -v "CONDA_BACKUP_" | sort
    rm -f /tmp/old-env-$$.txt /tmp/new-env-$$.txt || true
  fi
fi
