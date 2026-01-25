<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Gene Ontology -> Unified Gene Function Ontology Schema Mapping
  ============================================================
  Source: ./schema.json (Gene Ontology terms and annotations)
  Target: ../schema.json (Unified 4.5 Gene Function Ontology Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ goid                        │ term_id                 │ Direct       │
  │ label                       │ term_name               │ Direct       │
  │ namespace                   │ ontology_aspect         │ Direct       │
  │ aspect                      │ aspect_code             │ Direct       │
  │ definition                  │ definition              │ Direct       │
  │ synonyms                    │ synonyms                │ Direct       │
  │ is_obsolete                 │ is_obsolete             │ Direct       │
  │ relationships               │ relationships           │ Direct       │
  │ annotation                  │ go_annotation           │ Restructure  │
  │ gocam                       │ go_gocam                │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Missing optional fields -> null

  NOTES:
  - GO IDs follow GO:####### pattern
  - Three namespaces: biological_process, molecular_function, cellular_component
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

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- GO term identifier -->
      <fn:string key="term_id">
        <xsl:value-of select="fn:string[@key='goid']"/>
      </fn:string>

      <!-- Term name -->
      <fn:string key="term_name">
        <xsl:value-of select="fn:string[@key='label']"/>
      </fn:string>

      <!-- Namespace (aspect) -->
      <fn:string key="ontology_aspect">
        <xsl:value-of select="fn:string[@key='namespace']"/>
      </fn:string>

      <!-- Single-letter aspect code -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='aspect'] and fn:string[@key='aspect'] != ''">
          <fn:string key="aspect_code">
            <xsl:value-of select="fn:string[@key='aspect']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="aspect_code"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Definition -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='definition'] and fn:string[@key='definition'] != ''">
          <fn:string key="definition">
            <xsl:value-of select="fn:string[@key='definition']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="definition"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Is obsolete -->
      <xsl:choose>
        <xsl:when test="fn:boolean[@key='is_obsolete']">
          <fn:boolean key="is_obsolete">
            <xsl:value-of select="fn:boolean[@key='is_obsolete']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:boolean key="is_obsolete">false</fn:boolean>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Synonyms -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='synonyms']">
          <fn:array key="synonyms">
            <xsl:for-each select="fn:array[@key='synonyms']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synonyms"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Relationships -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='relationships']">
          <fn:map key="relationships">
            <xsl:if test="fn:map[@key='relationships']/fn:array[@key='is_a']">
              <fn:array key="is_a">
                <xsl:for-each select="fn:map[@key='relationships']/fn:array[@key='is_a']/fn:string">
                  <fn:string><xsl:value-of select="."/></fn:string>
                </xsl:for-each>
              </fn:array>
            </xsl:if>
            <xsl:if test="fn:map[@key='relationships']/fn:array[@key='part_of']">
              <fn:array key="part_of">
                <xsl:for-each select="fn:map[@key='relationships']/fn:array[@key='part_of']/fn:string">
                  <fn:string><xsl:value-of select="."/></fn:string>
                </xsl:for-each>
              </fn:array>
            </xsl:if>
            <xsl:if test="fn:map[@key='relationships']/fn:array[@key='regulates']">
              <fn:array key="regulates">
                <xsl:for-each select="fn:map[@key='relationships']/fn:array[@key='regulates']/fn:string">
                  <fn:string><xsl:value-of select="."/></fn:string>
                </xsl:for-each>
              </fn:array>
            </xsl:if>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="relationships"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GO-SPECIFIC FIELDS ========== -->

      <!-- GO annotation -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='annotation']">
          <fn:map key="go_annotation">
            <xsl:copy-of select="fn:map[@key='annotation']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="go_annotation"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- GO-CAM model -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='gocam']">
          <fn:map key="go_gocam">
            <xsl:copy-of select="fn:map[@key='gocam']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="go_gocam"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Gene Ontology</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='goid']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Convert empty/placeholder to empty string (for null handling) -->
  <xsl:function name="local:nullify" as="xs:string">
    <xsl:param name="value" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>
