import pytest
import subprocess
import os
 
import pysam

# Constants
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = 'merge_duplicates.py'
SCRIPT = os.path.join(TEST_DIR, '..', 'scripts', SCRIPT_NAME)
INPUT_DIR = os.path.join(TEST_DIR, 'merge_duplicates', 'input')
TRUTH_DIR = os.path.join(TEST_DIR, 'merge_duplicates', 'truth')
OUTPUT_DIR = os.path.join(TEST_DIR, 'merge_duplicates', 'output')
TYPE = ['position', 'mis_as_ref', 'repeat', 'none', 'warning']
ERROR_TYPE = ['right_shift_error']

# Run command
def run_script(vcf_type: str) -> None:

    invcf = os.path.join(INPUT_DIR, vcf_type + '.input.vcf.gz')
    outvcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')

    # create output dir if not exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    # clean previous test run if exist
    if os.path.exists(outvcf):
        os.remove(outvcf)

    command = [SCRIPT, '-i', invcf, '-o',  outvcf]
    
    if vcf_type == 'position':
        command += ['-c', 'position','--track', 'AT']
    elif vcf_type == 'mis_as_ref':
        command += ['-c', 'position','--track', 'AT', '--merge-mis-as-ref']
    elif vcf_type == 'repeat':
        command += ['-c', 'repeat', '--max-repeat', '100', '--track', 'ID', '--keep-order']
    elif vcf_type == 'none':
        command += ['-c', 'none', '--keep-order']
    elif vcf_type == 'warning':
        command += ['-c', 'position', '--keep-order']

    subprocess.run(command, check=True)


# Test if the command runs successfully
@pytest.mark.order(1)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_script_execution(vcf_type: str) -> None:
    try:
        run_script(vcf_type)  # This will raise CalledProcessError on failure
        assert True
    except subprocess.CalledProcessError as e:
        pytest.fail(f"{SCRIPT_NAME} failed with error: {e}")


# Compare output files with truth
@pytest.mark.order(2)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_output(vcf_type: str) -> None:

    # Get truth and test files for this suffix
    file_truth_vcf = os.path.join(TRUTH_DIR, vcf_type + '.output.vcf.gz')
    file_output_vcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')

    # Read files
    truth_vcf = pysam.VariantFile(file_truth_vcf, 'rb')
    output_vcf = pysam.VariantFile(file_output_vcf, 'rb')

    # Compare header
    for truth_header, output_header in zip(truth_vcf.header.records, output_vcf.header.records):
        assert str(truth_header) == str(output_header), f"Header mismatch for VCF: {vcf_type}"
    
    # Compare samples
    assert list(truth_vcf.header.samples) == list(output_vcf.header.samples), f"Sample mismatch for VCF: {vcf_type}"

    # Compare records
    for truth_record, output_record in zip(truth_vcf, output_vcf):
        assert str(truth_record) == str(output_record), f"Record mismatch for VCF: {vcf_type}"

    truth_vcf.close()
    output_vcf.close()

    os.remove(file_output_vcf)

# Catch error
@pytest.mark.order(3)
@pytest.mark.parametrize("error_type", ERROR_TYPE)
def test_error(error_type: str) -> None:
    
    invcf = os.path.join(INPUT_DIR, error_type + '.input.vcf.gz')
    outvcf = os.path.join(OUTPUT_DIR, error_type + '.output.vcf.gz')

    command = [SCRIPT, '-i', invcf, '-o',  outvcf]

    if error_type == 'right_shift_error':
        expected_error = 'ValueError: Cannot right shift'

    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        subprocess.run(command, check=True, capture_output=True, text=True)

    assert expected_error in excinfo.value.stderr

    os.remove(outvcf)

# A concatenated REF is synthesised by joining indel sequences, not read from the reference, so
# two concatenations at one position can each agree with every input record and still disagree
# with one another. `bcftools norm -m +any` rejects such a pair outright ("The REF prefixes
# differ") and takes the surrounding pipe down with it, so we must never emit one.
@pytest.mark.order(4)
def test_concat_refs_mutually_compatible(tmp_path) -> None:

    # the reference starts ACTCTC..., so every REF below is a genuine prefix of it, exactly as
    # `bcftools norm -f` would leave them. sample1 carries v1+v2, which concatenates to ACTCTCT;
    # sample2 carries v1+v3, which concatenates to ACTCTCCT. Each is compatible with all three
    # inputs, but they diverge from each other at offset 6.
    records = [('ACT', '1', '1'), ('ACTCT', '1', '0'), ('ACTCTC', '0', '1')]
    header = ['##fileformat=VCFv4.2',
              '##contig=<ID=c,length=1000>',
              '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
              '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
              '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">',
              '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2']

    invcf = tmp_path / 'ref_conflict.input.vcf'
    outvcf = tmp_path / 'ref_conflict.output.vcf.gz'
    with open(invcf, 'w') as in_file:
        in_file.write('\n'.join(header) + '\n')
        for i, (ref, gt1, gt2) in enumerate(records, 1):
            in_file.write('c\t10\tv{}\t{}\tA\t60\t.\tAC={};AF=0.5;AN=2;NS=2\tGT\t{}\t{}\n'.format(
                i, ref, int(gt1) + int(gt2), gt1, gt2))

    subprocess.run([SCRIPT, '-i', str(invcf), '-o', str(outvcf)], check=True)

    refs_by_pos = {}
    output_vcf = pysam.VariantFile(str(outvcf))
    for record in output_vcf:
        refs_by_pos.setdefault((record.chrom, record.pos), []).append(record.ref)
    output_vcf.close()

    for (chrom, pos), refs in refs_by_pos.items():
        for i, ref_a in enumerate(refs):
            for ref_b in refs[i + 1:]:
                shorter, longer = sorted((ref_a, ref_b), key=len)
                assert longer.startswith(shorter), \
                    f"incompatible REFs at {chrom}:{pos}: {ref_a} vs {ref_b}"
