# Unmodified function from T0hid/hpo-classification-agent at 03bc5f6be6456a3ca2c5206e399f0a38f879ce57
# MIT; see upstream_LICENSE.txt. Synthetic parity reference only.

def calculate_severity(row):
    """Severity classification rules based on trait counts per tier."""
    count_tier_1 = row.get('1', 0)
    count_tier_2 = row.get('2', 0)
    count_tier_3 = row.get('3', 0)
    if count_tier_1 > 1:  return 'Profound'
    if count_tier_1 == 1: return 'Severe'
    if count_tier_2 >= 1:
        return 'Severe' if (count_tier_2 + count_tier_3) >= 4 else 'Moderate'
    if count_tier_3 >= 1: return 'Moderate'
    return 'Mild'
