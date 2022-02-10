# This is a collection of scripts that will allow manipulation of CAMI profiling files
import numpy as np
from src.utils.ProfilingToolsABC import ProfileABC


# TODO: make sure that I'm not deleting the root "-1" that way Unifrac picks up on the missing superkingdoms

class Profile(ProfileABC):
    def __init__(self, sample_metadata=None, profile=None, branch_length_fun=lambda x: 1/x):
        super().__init__(sample_metadata, profile)
        self.root_len = 1  # the length you want between the "root" of "-1" and the superkingdom level (eg. Bacteria)
        self.branch_len_func = branch_length_fun  # Given a node n at depth d in the tree, branch_len_func(d)
        # is how long you want the branch length between n and ancestor(n) to be
        self._data["-1"]["branch_length"] = self.root_len
        self.parse_file()  # TODO: this sets all the branch lengths to 1 currently

    def parse_file(self):
        _data = self._data
        _all_keys = self._all_keys
        _header = self._header
        for k, v in self.sample_metadata.items():
            _header.append('{}:{}'.format(k, v))

        # populate all the correct keys
        for prediction in self.profile:
            _all_keys.append(prediction.taxid.strip())

        # crawl over all profiles tax_path and create the ancestors and descendants list
        for prediction in self.profile:
            tax_id = prediction.taxid.strip()
            tax_path = prediction.taxpath.strip().split("|")  # this will be a list, join up late
            if tax_id not in _data:
                _data[tax_id] = dict()
            else:
                raise Exception(f"Improperly formatted profile: row starting with {tax_id} shows up more than once")
            _data[tax_id]["tax_path"] = tax_path

            # populate abundance
            _data[tax_id]["abundance"] = prediction.percentage

            # populate tax path sn
            if not (prediction.taxpathsn is None):  # might not be present
                _data[tax_id]["tax_path_sn"] = prediction.taxpathsn.strip().split("|")  # this will be a list, join up later

            # populate the rank
            _data[tax_id]["rank"] = prediction.rank.strip()

            # populate the branch length
            _data[tax_id]["branch_length"] = self.tax_path_to_branch_len(tax_path, self.branch_len_func, self.root_len)

            # Find the ancestors
            if len(tax_path) <= 1:  # note, due to the format, we will never run into the case tax_path == []
                _data[tax_id]["ancestor"] = "-1"  # no ancestor, it's a root
            else:  # go from the bottom up, looking for an ancestor that is an acceptable key
                ancestor = "-1"  # this is the default
                tax_path_rev = tax_path[::-1]
                for potential_ancestor in tax_path_rev:
                    if potential_ancestor != tax_id and potential_ancestor in _all_keys:
                        ancestor = potential_ancestor
                        break  # you found the ancestor, so can quit looking
                _data[tax_id]["ancestor"] = ancestor

            # Create a placeholder descendant key initialized to [], just so each tax_id has a descendant key associated to it
            if "descendants" not in _data[tax_id]:  # if this tax_id doesn't have a descendant list,
                _data[tax_id]["descendants"] = list()  # initialize to empty list

        self._add_descendants()
        self._delete_missing()  # make sure there aren't any missing internal nodes

    def _add_descendants(self):
        """
        Idea here is to look at all the ancestors of each key, and make the key the descendant of that ancestor
        Returns
        -------
        None: modifies Profile in place
        """
        _data = self._data
        _all_keys = self._all_keys
        for prediction in self.profile:
            tax_id = prediction.taxid.strip()  # the tax ID we are looking at
            ancestor = _data[tax_id]['ancestor']  # the tax ID's ancestor
            if tax_id not in _data[ancestor]['descendants']:
                _data[ancestor]['descendants'].append(tax_id)  # so make the tax ID we're looking at the descendant of the ancestor

    def _delete_missing(self):
        """
        Deletes from the descendants all those taxids that aren't keys in the profile (i.e. there is no line that starts with that taxID)
        Returns
        -------
        none: modifies Profile in place
        """
        for key in self._data:
            clean_descendants = []
            for descendant in self._data[key]["descendants"]:
                if descendant in self._all_keys:  # if it's one of the taxids that the line starts with, add it
                    clean_descendants.append(descendant)
                else:
                    pass  # don't include the taxids that aren't actually in the final tax tree
            self._data[key]["descendants"] = clean_descendants
        return

    @staticmethod
    def tax_path_to_branch_len(tax_path, func, root_len=1):
        """
        This function modifies the branch lengths based on the input tax_path.
        intent is: ["2", "", "123", "456"] would result in a branch length of func(4)
        Parameters
        ----------
        tax_path : a list of strings (tax ID's)
        func : a function whose argument is the depth in the tree of a tax ID, and whose output is the branch length
               from the tax ID to its ancestor.
        root_len : how long you want the root of the tree "-1" to be to the descendants (eg. "-1" -> "Bacteria")
        Returns
        -------
        float
        """
        # eg. "-1" -> "Bacteria" should have a branch length of root_len
        if not tax_path:
            return root_len
        else:
            depth_in_tree = len(tax_path)  # this takes into account that the tax_path doesn't include the root of "-1"
            return func(depth_in_tree)

    def make_unifrac_input_and_normalize(self, other):
        if not isinstance(other, Profile):
            raise Exception
        _data = self._data
        _other_data = other._data

        _data_keys = _data.keys()
        tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
        _other_data_keys = _other_data.keys()
        tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
        tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
        all_keys = set(_data_keys)
        all_keys.update(_other_data_keys)  # all the taxID's in the union of self and other profile
        nodes_in_order = []
        for path_length in range(tax_path_lengths, 0, -1):
            for key in all_keys:
                if key in _data:
                    if len(_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
                elif key in _other_data:
                    if len(_other_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
        # Make the graph
        # Put the root at the very end
        if '-1' in nodes_in_order:
            nodes_in_order.pop(nodes_in_order.index('-1'))
            nodes_in_order.append('-1')
        else:
            nodes_in_order.append('-1')
        Tint = dict()
        lint = dict()
        for key in nodes_in_order:
            if key in _data:
                if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
                    ancestor = _data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _data[key]["branch_length"]
            elif key in _other_data:
                if "ancestor" in _other_data[key]:
                    ancestor = _other_data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _other_data[key]["branch_length"]
        nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))  # maps '45202.15' -> 0 (i.e taxID to integer index)

        # Now need to change over to the integer-based indexing
        Tint2 = dict()
        lint2 = dict()
        nodes_in_order2 = []
        for key in nodes_in_order:
            if key in Tint:
                ancestor = Tint[key]
                Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
                if (key, ancestor) in lint:
                    lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
            nodes_in_order2.append(nodes_to_index[key])

        # Next make the probability distributions
        # Would be nice if I could find a non-destructive way to subtract up and normalize

        # Do it for P
        self._subtract_down()
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
        P = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _data:
                P[key_ind] = _data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] *= 100
        self._add_up()

        # Next do for Q
        other._subtract_down()
        keys = _other_data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _other_data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
        Q = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _other_data:
                Q[key_ind] = _other_data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        for key in keys:
            if total_abundance > 0:
                _other_data[key]["abundance"] *= 100
        other._add_up()

        return Tint2, lint2, nodes_in_order2, P, Q

    def make_unifrac_input_no_normalize(self, other):
        if not isinstance(other, Profile):
            raise Exception
        _data = self._data
        _other_data = other._data

        _data_keys = _data.keys()
        tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
        _other_data_keys = _other_data.keys()
        tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
        tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
        all_keys = set(_data_keys)
        all_keys.update(_other_data_keys)  # all the taxID's in the union of self and other profile
        nodes_in_order = []
        for path_length in range(tax_path_lengths, 0, -1):
            for key in all_keys:
                if key in _data:
                    if len(_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
                elif key in _other_data:
                    if len(_other_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
        # Make the graph
        # Put the root at the very end
        if '-1' in nodes_in_order:
            nodes_in_order.pop(nodes_in_order.index('-1'))
            nodes_in_order.append('-1')
        else:
            nodes_in_order.append('-1')
        Tint = dict()
        lint = dict()
        for key in nodes_in_order:
            if key in _data:
                if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
                    ancestor = _data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _data[key]["branch_length"]
            elif key in _other_data:
                if "ancestor" in _other_data[key]:
                    ancestor = _other_data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _other_data[key]["branch_length"]
        nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))  # maps '45202.15' -> 0 (i.e taxID to integer index)

        # Now need to change over to the integer-based indexing
        Tint2 = dict()
        lint2 = dict()
        nodes_in_order2 = []
        for key in nodes_in_order:
            if key in Tint:
                ancestor = Tint[key]
                Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
                if (key, ancestor) in lint:
                    lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
            nodes_in_order2.append(nodes_to_index[key])

        # Next make the probability distributions
        # Would be nice if I could find a non-destructive way to subtract up and normalize

        # Do it for P
        self._subtract_down()
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                #_data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
                pass
        P = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _data:
                P[key_ind] = _data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        #for key in keys:
        #    if total_abundance > 0:
        #        _data[key]["abundance"] *= 100
        self._add_up()

        # Next do for Q
        other._subtract_down()
        keys = _other_data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _other_data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                #_other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
                pass
        Q = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _other_data:
                Q[key_ind] = _other_data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        #for key in keys:
        #    if total_abundance > 0:
        #        _other_data[key]["abundance"] *= 100
        other._add_up()

        return Tint2, lint2, nodes_in_order2, P/100., Q/100.


def test_normalize():
    import EMDUnifrac as EMDU
    from load_data import open_profile_from_tsv
    import os
    # test files
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/agitated_blackwell_7"
    file_path2 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/goldstandard_low_1.bin"

    # import one test profile
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile1 = profile_list[0]
    profile1 = Profile(sample_metadata=metadata1, profile=profile1)

    # import another test profile
    profile_list = open_profile_from_tsv(file_path2, False)
    name2, metadata2, profile2 = profile_list[0]
    profile2 = Profile(sample_metadata=metadata2, profile=profile2)

    print("Normalized:")
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile2)  # normalized
    print(f"P sum: {np.sum(P)}")
    print(f"Q sum: {np.sum(Q)}")
    (weighted_norm, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    (unweighted_norm, diffab) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    print(f"weighted, normalized: {weighted_norm}")
    print(f"unweighted, normalized: {unweighted_norm}")


    print("No normalized:")
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile1 = profile_list[0]
    profile1 = Profile(sample_metadata=metadata1, profile=profile1)
    # import another test profile
    profile_list = open_profile_from_tsv(file_path2, False)
    name2, metadata2, profile2 = profile_list[0]
    profile2 = Profile(sample_metadata=metadata2, profile=profile2)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_no_normalize(profile2)  # not normalized
    print(f"P sum: {np.sum(P)}")
    print(f"Q sum: {np.sum(Q)}")
    P_missing_mass = 1-np.sum(P)
    (weighted_no_norm, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    (unweighted_no_norm, diffab) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    print(f"weighted, not normalized: {weighted_no_norm}")
    print(f"unweighted, not normalized: {unweighted_no_norm}")
    print(f"weighted, not normalized, but missing mass added to root node: {weighted_no_norm + P_missing_mass}")

    assert unweighted_norm == unweighted_no_norm
    assert weighted_norm != unweighted_no_norm
    return


def test_branch_lengths_all_1():
    from load_data import open_profile_from_tsv
    import os
    # test file
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/small1.profile"
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile_fernando = profile_list[0]

    # Test with branch lengths of 1
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: 1)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    #print(f"Tint_new: {Tint_new}")
    #print(f"lint_new: {lint_new}")
    #print(f"P: {P}")
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    for val in lint.values():
        assert val == 1
    correct_vals = {'0': 0.20, '1': 0.50, '2': 0.20, '3': 0.0, '4': 0.10, '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val

    # test with branch lengths of x
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: x)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    correct_lints = {('1', '3'): 3, ('3', '4'): 2, ('4', '-1'): 1, ('2', '4'): 3, ('0', '5'): 5, ('5', '-1'): 1}
    for key, val in correct_lints.items():
        assert lint_new[key] == correct_lints[key]
    correct_vals = {'0': 0.20, '1': 0.50, '2': 0.20, '3': 0.0, '4': 0.10, '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val


def test_branch_lengths_all_2():
    from load_data import open_profile_from_tsv
    import os
    # test file
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/small2.profile"
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile_fernando = profile_list[0]

    # Test with branch lengths of 1
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: 1)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    #print(f"Tint_new: {Tint_new}")
    #print(f"lint_new: {lint_new}")
    #print(f"P: {P}")
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    for val in lint.values():
        assert val == 1
    correct_vals = {'0': 1/9., '1': 5/9., '2': 2/9., '3': 0.0, '4': 1/9., '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val

    # test with branch lengths of x
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: x)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    correct_lints = {('1', '3'): 3, ('3', '4'): 2, ('4', '-1'): 1, ('2', '4'): 3, ('0', '5'): 5, ('5', '-1'): 1}
    for key, val in correct_lints.items():
        assert lint_new[key] == correct_lints[key]
    correct_vals = {'0': 1/9., '1': 5/9., '2': 2/9., '3': 0.0, '4': 1/9., '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val


def test_no_normalize():
    from load_data import open_profile_from_tsv
    import os
    # test file
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/small2.profile"
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile_fernando = profile_list[0]

    # Test with branch lengths of 1
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: 1)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_no_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    #print(f"Tint_new: {Tint_new}")
    #print(f"lint_new: {lint_new}")
    #print(f"P: {P}")
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    for val in lint.values():
        assert val == 1
    correct_vals = {'0': 0.10, '1': 0.50, '2': 0.20, '3': 0.0, '4': 0.10, '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val


if __name__ == "__main__":
    test_normalize()
    test_branch_lengths_all_1()
    test_branch_lengths_all_2()
    test_no_normalize()
