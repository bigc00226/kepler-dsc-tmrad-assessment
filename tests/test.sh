#!/bin/bash
# Writes /logs/verifier/reward.txt (exactly 0 or 1) on every path.
mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt

cd /tests
python -m pytest test_outputs.py -rA -p no:cacheprovider \
    --ctrf /logs/verifier/ctrf.json > /logs/verifier/pytest.log 2>&1
rc=$?
cat /logs/verifier/pytest.log

if [ "$rc" -eq 0 ]; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
exit 0
