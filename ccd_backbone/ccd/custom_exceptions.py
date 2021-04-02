class DNASearchTimeoutError(Exception):
    pass

class IdNotFoundError(ValueError):
    pass

class SequenceNotFoundError(Exception):
    pass

class NotAnORF(Exception):
    pass

class IGiveUpError(Exception):
    pass

class NoSequenceFound(Exception):
    pass

class PicrNotResponding(Exception):
    pass

class NoIsoformsPresentError(Exception):
    pass

class NoNeedToAlignError(Exception):
    pass

class NoAlignmentCached(Exception):
    pass
