#!/usr/bin/env bash
# Complete workflow demonstration for wombat CLI

echo "=========================================="
echo "Wombat CLI - Complete Workflow Demo"
echo "=========================================="
echo ""

# 1. Show help
echo "1. Showing CLI help:"
echo "-------------------"
uv run wombat --help
echo ""

echo "2. Showing format command help:"
echo "------------------------------"
uv run wombat format --help
echo ""

# 2. Create test data if needed
if [ ! -f "tests/test_small.tsv" ]; then
    echo "3. Creating test data..."
    cat > tests/test_small.tsv << 'EOF'
CHROM	POS	REF	ALT	(null)	Sample1:GT	Sample2:GT	Sample3:GT
chr1	100	A	T	DP=30;AF=0.5;AC=2	0/1	1/1	0/0
chr1	200	G	C	DP=45;AF=0.25;AC=1	0/0	0/1	0/0
chr2	150	C	G	DP=60;AF=0.75;AC=3	1/1	0/1	1/1
EOF
fi

# 3. Format the test file
echo "3. Formatting test file with verbose output:"
echo "-------------------------------------------"
uv run wombat format tests/test_small.tsv -o tests/demo_output.tsv --verbose
echo ""

# 4. Show the output
echo "4. First 10 lines of output:"
echo "---------------------------"
head -10 tests/demo_output.tsv | column -t -s $'\t'
echo ""

# 5. Count rows
echo "5. Row counts:"
echo "-------------"
echo "Input rows:  $(tail -n +2 tests/test_small.tsv | wc -l)"
echo "Output rows: $(tail -n +2 tests/demo_output.tsv | wc -l)"
echo "(Output should have 3x input rows for 3 samples)"
echo ""

echo "=========================================="
echo "Demo complete!"
echo "=========================================="
