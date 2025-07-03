"""Helper functions to maintain consistent entity notation."""


def var_str(contig, pos, ref, alt):
    """Compose the variant ID string contig:pos:ref>alt"""
    return f'{contig}:{pos}:{ref}>{alt}'
