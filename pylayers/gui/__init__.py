"""
pylayers
=========

This file is adapted from scikit-learn package 

"""
import sys
__version__ = '0.12-git'

try:
    # This variable is injected in the __builtins__ by the build
    # process. It used to enable importing subpackages of sklearn when
    # the binaries are not built
    __PYLAYERS_SETUP__
except NameError:
    __PYLAYERS_SETUP__ = False

if __PYLAYERS_SETUP__:
    sys.stderr.write('Partial import of pylayers during the build process.\n')
    # We are not importing the rest of the scikit during the build
    # process, as it may not be compiled yet
else:
    try:
        from numpy.testing import nosetester

        class _NoseTester(nosetester.NoseTester):
            """ Subclass numpy's NoseTester to add doctests by default
            """

            def test(self, label='fast', verbose=1, extra_argv=['--exe'],
                            doctests=True, coverage=False):
                """Run the full test suite

                Examples
                --------
                This will run the test suite and stop at the first failing
                example
                >>> from pylayers import test
                >>> test(extra_argv=['--exe', '-sx']) #doctest: +SKIP
                """
                return super(_NoseTester, self).test(label=label, verbose=verbose,
                                        extra_argv=extra_argv,
                                        doctests=doctests, coverage=coverage)

        try:
            test = _NoseTester(raise_warnings="release").test
        except TypeError:
            # Older versions of numpy do not have a raise_warnings argument
            test = _NoseTester().test
        del nosetester
    except:
        pass

    __all__ = ['gis', 'signal', 'antprop', 'simul','util']

