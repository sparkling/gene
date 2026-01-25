<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ClassyFire -> Chemical Ontology Unified Schema Mapping
  ============================================================
  Source: ./schema.md (ClassyFire Chemical Classification)
  Target: ../schema.json (2.6 Chemical Ontology and Classification Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ chemont_id                  │ chemont_id              │ Direct       │
  │ smiles                      │ canonical_smiles        │ Direct       │
  │ inchikey                    │ inchi_key               │ Direct       │
  │ kingdom                     │ cf_kingdom              │ Object       │
  │ superclass                  │ cf_superclass           │ Object       │
  │ class                       │ cf_class                │ Object       │
  │ subclass                    │ cf_subclass             │ Object       │
  │ direct_parent               │ cf_direct_parent        │ Object       │
  │ alternative_parents         │ cf_alternative_parents  │ Object array │
  │ substituents                │ cf_substituents         │ String array │
  │ molecular_framework         │ cf_molecular_framework  │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - ClassyFire provides hierarchical chemical classification
  - Uses ChemOnt IDs (CHEMONTID:0000XXX format)
  - Kingdom -> Superclass -> Class -> Subclass -> Direct Parent
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
        <xsl:when test="local:is-null(fn:string[@key='chemont_id'])">
          <fn:null key="chemont_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="chemont_id">
            <xsl:value-of select="fn:string[@key='chemont_id']"/>
          </fn:string>
          <fn:string key="compound_id">
            <xsl:value-of select="fn:string[@key='chemont_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchikey'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchikey']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLASSYFIRE HIERARCHY ========== -->

      <xsl:choose>
        <xsl:when test="fn:map[@key='kingdom']">
          <xsl:copy-of select="fn:map[@key='kingdom']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_kingdom"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='superclass']">
          <xsl:copy-of select="fn:map[@key='superclass']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_superclass"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='class']">
          <xsl:copy-of select="fn:map[@key='class']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_class"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='subclass']">
          <xsl:copy-of select="fn:map[@key='subclass']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_subclass"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='direct_parent']">
          <xsl:copy-of select="fn:map[@key='direct_parent']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_direct_parent"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ALTERNATIVE PARENTS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='alternative_parents']/fn:map">
          <fn:array key="cf_alternative_parents">
            <xsl:for-each select="fn:array[@key='alternative_parents']/fn:map">
              <xsl:copy-of select="."/>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_alternative_parents"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SUBSTITUENTS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='substituents']/fn:string">
          <fn:array key="cf_substituents">
            <xsl:for-each select="fn:array[@key='substituents']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cf_substituents"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MOLECULAR FRAMEWORK ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molecular_framework'])">
          <fn:null key="cf_molecular_framework"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cf_molecular_framework">
            <xsl:value-of select="fn:string[@key='molecular_framework']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">ClassyFire</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='inchikey']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
