abstract type AbstractNode end

mutable struct Directive <: AbstractNode
    record::GFF3.Record
end

AbstractTrees.children(::Directive) = ()

mutable struct Comment <: AbstractNode
    record::GFF3.Record
end

AbstractTrees.children(::Comment) = ()

mutable struct Feature <: AbstractNode
    record::GFF3.Record
    seqid::Union{Missing, String}
    featuretype::Union{Missing, String}
    source::Union{Missing, String}
    seqstart::Union{Missing, Int}
    seqend::Union{Missing, Int}
    score::Union{Missing, Float64}
    strand::Union{Missing, GFF3.GenomicFeatures.Strand}
    phase::Union{Missing, Int}
    attributes::Vector{Pair{String,Vector{String}}}
    children::Vector{Feature}
end

GFF3.Record(f::Feature) = f.record
GFF3.seqid(f::Feature) = f.seqid
GFF3.featuretype(f::Feature) = f.featuretype
GFF3.source(f::Feature) = f.source
GFF3.seqstart(f::Feature) = f.seqstart
GFF3.seqend(f::Feature) = f.seqend
GFF3.score(f::Feature) = f.score
GFF3.strand(f::Feature) = f.strand
GFF3.phase(f::Feature) = f.phase
GFF3.attributes(f::Feature) = f.attributes 
AbstractTrees.children(f::Feature) = f.children

mutable struct Chromosome <: AbstractNode
    seqid::String
    children::Vector{Feature}
    allnode::Dict{String, Feature}
end

GFF3.seqid(c::Chromosome) = c.seqid
AbstractTrees.children(c::Chromosome) = c.children

AbstractTrees.nodevalue(n::AbstractNode) = n.record

AbstractTrees.printnode(io::IO, d::Directive) = print(io, d.record)
AbstractTrees.printnode(io::IO, c::Comment) = print(io, c.record)
AbstractTrees.printnode(io::IO, f::Feature) = print(io, f.record)
AbstractTrees.printnode(io::IO, c::Chromosome) = print(io, c.seqid)

Base.show(io::IO, c::Chromosome) = AbstractTrees.print_tree(io, c)
Base.show(io::IO, f::Feature) = print(io, f.record)

function Chromosome(c::String)
    children = Vector{Feature}[]
    allnodes = Dict{String,Feature}()
    return Chromosome(c, children, allnodes)
end

function subscribe!(c::Chromosome, f::Feature)
    if hasid(f)
        ids = getids(f)
        for id in ids 
            push!(c.allnode, Pair(id, f))
        end
    end
end

getfeature(c::Chromosome, id::AbstractString, m) = get(c.allnode, id, m)

function Feature(record::GFF3.Record)
    seqid = GFF3.hasseqid(record) ? GFF3.seqid(record) : missing
    featuretype = GFF3.hasfeaturetype(record) ? GFF3.featuretype(record) : missing
    source = GFF3.hassource(record) ? GFF3.source(record) : missing
    seqstart = GFF3.hasseqstart(record) ? GFF3.seqstart(record) : missing
    seqend = GFF3.hasseqend(record) ? GFF3.seqend(record) : missing
    score = GFF3.hasscore(record) ? GFF3.score(record) : missing
    strand = GFF3.hasstrand(record) ? GFF3.strand(record) : missing
    phase = GFF3.hasphase(record) ? GFF3.phase(record) : missing
    attributes = GFF3.attributes(record)
    children = Vector{Feature}[]
    return Feature(record, seqid, featuretype, source, seqstart, seqend, score, strand,  phase, attributes, children)
end

function update!(f::Feature)
    attributes_str = join(map(attr -> first(attr) * "=" * join(map(x -> URIParser.escape(x), last(attr)), ","), f.attributes), ";")
    seqid = encode(f.seqid)
    featuretype = encode(f.featuretype)
    source = encode(f.source)
    seqstart = encode(f.seqstart)
    seqend = encode(f.seqend)
    score = encode(f.score)
    strand = encode(f.strand)
    phase = encode(f.phase)
    attributes = encode(attributes_str)

    content = join(
        [seqid,
        source,
        featuretype,
        seqstart,
        seqend,
        score,
        strand,
        phase,
        attributes
        ],
        "\t"
    )
    new_record = GFF3.Record(content)
    f.record = new_record
end

function add_child!(c::Chromosome, child::Feature)
    push!(c.children,child)
end

function add_child!(f::Feature, child::Feature)
    push!(f.children,child)
end

encode(m::Missing) = "."
encode(s::String) = isempty(s) ? "." : s
encode(n::Number) = string(n) 
function encode(g::GFF3.GenomicFeatures.Strand)
    if g == GFF3.GenomicFeatures.STRAND_NA
        return "?"
    elseif g == GFF3.GenomicFeatures.STRAND_POS
        return "+"
    elseif g == GFF3.GenomicFeatures.STRAND_NEG
        return "-"
    elseif g ==GFF3.GenomicFeatures.STRAND_BOTH
        return "."
    end
end

function getids(f::Feature)
    return get(Dict(f.attributes), "ID", AbstractString[])
end

function getparentids(f::Feature)
    return get(Dict(f.attributes), "Parent", AbstractString[])
end

hasid(f::Feature) = !isempty(getids(f))
hasparent(f::Feature) = !isempty(getparentids(f))

function replaceid!(f::Feature, old_new::Pair{String,String})
    previous_ids = getids(f)
    if isempty(previous_ids)
        @warn "No id attributes"
        return
    end
    
    if !in(first(old_new), previous_ids)
        @warn "No old id in previous attributes"
        return
    end
    new_ids = replace(previous_ids, old_new)
    new_attributes = replace(f.attributes, Pair("ID", previous_ids)=>Pair("ID", new_ids))
    f.attributes = new_attributes

    for c in children(f)
        previous_parents = getparentids(c)
        if isempty(previous_parents)
            @warn "No parent attributes"
            continue
        end
    
        if !in(first(old_new), previous_parents)
            @warn "No old parent in previous attributes"
        end
        new_parents = replace(previous_parents, old_new)
        new_attributes = replace(c.attributes, Pair("Parent", previous_parents)=>Pair("Parent", new_parents))
        c.attributes = new_attributes
    end
end

function getnames(f::Feature)
    return get(Dict(f.attributes), "Name", AbstractString[])
end

hasname(f::Feature) = !isempty(getnames(f))

function replacename!(f::Feature, old_new::Pair{String, String})
    previous_names = getnames(f)
    if isempty(previous_names)
        @warn "No name attributes"
        return
    end

    if !in(first(old_new), previous_names)
        @warn "No old name in previous attributes"
    end

    new_names = replace(previous_names, old_new)
    new_attributes = replace(f.attributes, Pair("Name", previous_names)=>Pair("Name", new_names))
    f.attributes = new_attributes
end

function parse(reader::GFF3.Reader)
    chromosomes = Dict{String, Chromosome}()
    chr = Chromosome("")
    for record in reader
        if GFF3.isdirective(record)
            d = Directive(record)
            #nothing / TODO do something
            continue
        elseif GFF3.iscomment(record)
            c = Comment(record)
            #nothing / TODO do something
            continue
        elseif GFF3.isfeature(record)
            #nothing happen
        else
            #not reach
        end

        feature = Feature(record)
    
        if chr.seqid != feature.seqid
            chr = Chromosome(GFF3.seqid(record))
            push!(chromosomes, GFF3.seqid(chr)=>chr)
        end

        subscribe!(chr, feature)

        parents = getparentids(feature)
        if isempty(parents)
            add_child!(chr, feature)
            continue
        end
    
        for p in parents
            parent_feature = getfeature(chr, p, missing)
            if ismissing(parent_feature)
                @warn "parent feature does not exist"
                continue
            else
                add_child!(parent_feature, feature)
            end
        end
    end
    return chromosomes
end

function updateall!(c::Chromosome)
    for child in children(c)
        updateall!(child)
    end
end

function updateall!(f::Feature)
    update!(f)
    for child in children(f)
        updateall!(child)
    end
end

function Base.write(writer::GFF3.Writer, c::Chromosome)
    for child in AbstractTrees.children(c)
        write(writer, child)
    end 
end

function Base.write(writer::GFF3.Writer, f::Feature)
    write(writer, f.record)
    for child in AbstractTrees.children(f)
        write(writer, child)
    end 
end