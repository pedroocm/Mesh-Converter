import color_map

def test_get_tag():
    assert color_map.get_tag_number(200) == 2
    assert color_map.get_tag_number(4500.23) == 3
    assert color_map.get_tag_number(95320.1) == 4
    assert color_map.get_tag_number(754321.245) == 5
    assert color_map.get_tag_number(3.4567) == 0
    assert color_map.get_tag_number(11.50192) == 1

def test_get_tag_array():
    assert color_map.get_tag_array([200, 4500.23, 95320.1, 754321.245, 3.4567, 15.50192]) == [2, 3, 4, 5, 0, 1]