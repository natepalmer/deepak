
class Filters(dict):
    """ Dict-like object specifying the filters associated with an Experiment.

    Each key represents a filter name which must be the name of an attribute of a pafparser.PafRecord object.
    Each value must be a length 2 tuple where:
      the first member is a value to compare to the PafRecord attribute,
      the second member is the method by which to compare it, must be either:
        a string specifying one of the built-in comparison methods for type(value), e.g. 'eq' or 'lt',
        a function which takes the value and the attribute as arguments and returns the True/False result of the test
    """

    def __init__(self, *args, **kwargs):
        super(Filters, self).__init__(*args, **kwargs)

    def test(self, record, tests="all"):
        """ Test if record passes filters. Returns a length 2 tuple where the first member is True/False whether the
        record passes, and the second member is the name of the failed filter if applicable, or None if the record
        passes all filters.

        *record* is a PafRecord object to be tested
        """
        if tests == "all":
            tests = self.keys()
        elif type(tests) not in (list, tuple, set):
            tests = [tests]

        for k in tests:
            v = self[k]
            subject = getattr(record, k)
            if callable(v[1]):
                if not v[1](v[0], subject):
                    return False, k
            else:
                func = getattr(subject, "__"+v[1]+"__")
                if not func(v[0]):
                    return False, k

        return True, None
