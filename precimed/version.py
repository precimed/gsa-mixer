_MAJOR = "2"
_MINOR = "2"
# On main and in a nightly release the patch should be one ahead of the last
# released build.
_PATCH = "0"
# This is mainly for nightly builds which have the suffix ".dev$DATE". See
# https://semver.org/#is-v123-a-semantic-version for the semantics.
#_SUFFIX = ""
_SUFFIX = "dev-2025-09-09"

VERSION_SHORT = "{0}.{1}".format(_MAJOR, _MINOR)
VERSION = "{0}.{1}.{2}{3}".format(_MAJOR, _MINOR, _PATCH, _SUFFIX)

MASTHEAD = "***********************************************************************\n"
MASTHEAD += f"* (c) 2016-2025 MiXeR v{VERSION} software\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* Center for Multimodal Imaging and Genetics / UCSD\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

# TBD: include git commit hash via $(git rev-parse HEAD)
