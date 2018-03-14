import xmltree, os, sequtils, strutils, strtabs, xmlparser

type
    Allele = object
        name: string
        gGroup: string
        pGroup: string
        nucSequence: string

const outputTmpl = "$#\t$#\t$#\t$#\t$#\t$#"

template log(s: string) =
    when false:
        stderr.writeLine(s)
        stderr.flushFile

proc initAllele(data: XmlNode): Allele =
    ## Initialize allele
    result.name = data.attrs["name"]
    for allelenode in data:
        case allelenode.tag
        of "hla_g_group":
            result.gGroup = allelenode.attrs["status"]
        of "hla_p_group":
            result.pGroup = allelenode.attrs["status"]
        of "sequence":
            for feature in allelenode:
                if feature.tag == "nucsequence":
                    result.nucSequence = feature.innerText
        else:
            discard


proc getAlleles(hla: XmlNode, startpattern: string): seq[Allele] =
    result = @[]
    for node in hla:
        let name = node.attrs["name"]
        if name.startsWith(startpattern):
            let allele = initAllele(node)
            result.add allele

proc printHeader() =
    ## Print header
    echo(outputTmpl % [
        "Name", "fwd_seq", "rev_seq", "length", "Ggroup", "Pgroup"
    ])


proc printPrimerMatches(alleles: seq[Allele], fwd, rev: string) =
    ## Print alleles that bind primers
    for allele in alleles:
        let fwdindex = allele.nucSequence.find(fwd)
        if fwdindex == -1:
            continue
        let revindex = allele.nucSequence.find(rev)
        if revindex == -1:
            continue
        var length = abs(revindex - fwdindex)
        if fwdindex < revindex:
            length += rev.len
        else: length += fwd.len
        echo(outputTmpl % [allele.name,
                           allele.nucSequence[fwdindex ..< fwdindex+fwd.len],
                           allele.nucSequence[revindex ..< revindex+rev.len],
                           $length,
                           allele.gGroup,
                           allele.pGroup])


proc main() =
    const usage = """
usage: hlapcr hla.xml prefix fwd rev

Virtual PCR for HLA alleles.

Positional arguments:
    hla.xml         Path to hla.xml
    prefix          Prefix for alleles (ex. HLA-DQB1)
    fwd             Forward primer (5'-3')
    rev             Reverse primer (5'-3')

Primer examples:
    DQB1*02 F       CGTGCGTCTTGTGAGCAGAA
    DQB1-234R *VIC  AGTTGGAGCTCCGCACGAC

    DQB1*0602 F     CCCGCAGAGGATTTCGTGTT
    DQB1-56R *VIC   AGGAGAGGTGAGCGTCGTCG
"""
    if paramCount() < 4:
        quit(usage)

    let
        prefix = paramStr(2)
        fwd = paramStr(3).toUpperAscii()
        rev = paramStr(4).toUpperAscii()
    log("prefix: " & prefix)
    log("fwd primer: " & fwd)
    log("rev primer: " & rev)

    log "parsing hla xml file"
    let hla = loadXml(paramStr(1))

    log "locating alleles"
    let subset = getAlleles(hla, prefix)

    if subset.len == 0:
        quit("Found no alleles matching the prefix")

    printHeader()
    printPrimerMatches(subset, fwd, rev)

main()