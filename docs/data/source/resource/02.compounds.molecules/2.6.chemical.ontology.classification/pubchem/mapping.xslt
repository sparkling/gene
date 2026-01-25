<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PubChem -> Chemical Ontology Unified Schema Mapping
  ============================================================
  Source: ./schema.md (PubChem Chemical Database)
  Target: ../schema.json (2.6 Chemical Ontology and Classification Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ CID                         │ pubchem_cid             │ Direct       │
  │ Title                       │ name                    │ Direct       │
  │ CanonicalSMILES             │ canonical_smiles        │ Direct       │
  │ IsomericSMILES              │ isomeric_smiles         │ Direct       │
  │ InChI                       │ inchi                   │ Direct       │
  │ InChIKey                    │ inchi_key               │ Direct       │
  │ MolecularFormula            │ molecular_formula       │ Direct       │
  │ IUPACName                   │ iupac_name              │ Direct       │
  │ XLogP                       │ xlogp                   │ Direct       │
  │ ExactMass                   │ exact_mass              │ Direct       │
  │ TPSA                        │ tpsa                    │ Direct       │
  │ Complexity                  │ complexity              │ Direct       │
  │ bioactivity                 │ bioactivity             │ Object array │
  │ patents                     │ patents                 │ Object array │
  │ literature                  │ literature              │ Object array │
  │ Fingerprint2D               │ fingerprint2d           │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - PubChem uses CID (Compound ID) as primary identifier
  - Contains bioassay results, patents, and literature
  - 2D fingerprints for similarity searching
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='CID']">
          <fn:number key="pubchem_cid">
            <xsl:value-of select="fn:number[@key='CID']"/>
          </fn:number>
          <fn:string key="compound_id">
            <xsl:value-of select="fn:number[@key='CID']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pubchem_cid"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Title'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='Title']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='IUPACName'])">
          <fn:null key="iupac_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="iupac_name">
            <xsl:value-of select="fn:string[@key='IUPACName']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='CanonicalSMILES'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='CanonicalSMILES']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='IsomericSMILES'])">
          <fn:null key="isomeric_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="isomeric_smiles">
            <xsl:value-of select="fn:string[@key='IsomericSMILES']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='InChI'])">
          <fn:null key="inchi"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi">
            <xsl:value-of select="fn:string[@key='InChI']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='InChIKey'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='InChIKey']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MolecularFormula'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='MolecularFormula']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== COMPUTED PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='XLogP']">
          <fn:number key="xlogp">
            <xsl:value-of select="fn:number[@key='XLogP']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="xlogp"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='ExactMass']">
          <fn:number key="exact_mass">
            <xsl:value-of select="fn:number[@key='ExactMass']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="exact_mass"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='TPSA']">
          <fn:number key="tpsa">
            <xsl:value-of select="fn:number[@key='TPSA']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tpsa"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Complexity']">
          <fn:number key="complexity">
            <xsl:value-of select="fn:number[@key='Complexity']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="complexity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== FINGERPRINT ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Fingerprint2D'])">
          <fn:null key="fingerprint2d"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="fingerprint2d">
            <xsl:value-of select="fn:string[@key='Fingerprint2D']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BIOACTIVITY ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='bioactivity']/fn:map">
          <fn:array key="bioactivity">
            <xsl:for-each select="fn:array[@key='bioactivity']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='assay_id']">
                  <fn:string key="assay_id">
                    <xsl:value-of select="fn:string[@key='assay_id']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='activity']">
                  <fn:string key="activity">
                    <xsl:value-of select="fn:string[@key='activity']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='target']">
                  <fn:string key="target">
                    <xsl:value-of select="fn:string[@key='target']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="bioactivity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PATENTS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='patents']/fn:map">
          <fn:array key="patents">
            <xsl:for-each select="fn:array[@key='patents']/fn:map">
              <fn:map>
                <fn:string key="patent_id">
                  <xsl:value-of select="fn:string[@key='patent_id']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='title']">
                  <fn:string key="title">
                    <xsl:value-of select="fn:string[@key='title']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="patents"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== LITERATURE ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='literature']/fn:map">
          <fn:array key="literature">
            <xsl:for-each select="fn:array[@key='literature']/fn:map">
              <fn:map>
                <xsl:if test="fn:number[@key='pmid']">
                  <fn:number key="pmid">
                    <xsl:value-of select="fn:number[@key='pmid']"/>
                  </fn:number>
                </xsl:if>
                <xsl:if test="fn:string[@key='title']">
                  <fn:string key="title">
                    <xsl:value-of select="fn:string[@key='title']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="literature"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">PubChem</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='CID']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
