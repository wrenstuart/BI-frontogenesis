Was getting scalar index error when it tried to *checkpoint* ζ because it's an AbstractField, which is a subtype of AbstractArray, which it then tried to iterate over. This is a problem in the checkpointer only: the JLD2 writer can save ζ just fine.

Gemini says that the checkpointer saves ζ because it's an auxiliary field.

Question going forward for improving performance: can I save all the diagnostics I want just in tracked_fields rather than in auxiliary fields? If so, then I should be able to keep them in there as abstract fields (rather than computing them with Field()) and then I can save on a bunch of memory (and possibly boost performance too?)