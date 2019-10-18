function blocks = block_split(linestring)

pattern = '\*+ \<FLIGHT CONDITIONS AND REFERENCE QUANTITIES \*+';
blocks = regexp(linestring,pattern,'split');
blocks = blocks(2:end);