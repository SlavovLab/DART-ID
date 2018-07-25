<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
		xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:loc="http://regis-web.systemsbiology.net/pepXML">
  <xsl:output method="text" />

  <xsl:template match="/*/">
    <xsl:apply-templates select="loc:spectrum_query/loc:search_result/loc:search_hit" />
  </xsl:template>

  <xsl:template match="loc:spectrum_query/loc:search_result/loc:search_hit">
    <xsl:text>e </xsl:text>
    <xsl:value-of select="@peptide" /><xsl:text>&#10;</xsl:text>
    <xsl:text>c </xsl:text>
    <xsl:value-of select="../../@assumed_charge" /><xsl:text>&#10;</xsl:text>
    <xsl:text>r </xsl:text>
    <xsl:value-of select="@protein" /><xsl:text>&#10;</xsl:text>
    <xsl:for-each select="loc:alternative_protein">
      <xsl:text>r </xsl:text>
      <xsl:value-of select="@protein" /><xsl:text>&#10;</xsl:text>
    </xsl:for-each>
    <xsl:apply-templates select="loc:analysis_result/loc:peptideprophet_result" />
  </xsl:template>

<xsl:template match="loc:analysis_result/loc:peptideprophet_result">
  <xsl:text>p </xsl:text>
  <xsl:value-of select="@probability" /><xsl:text>&#10;</xsl:text>
</xsl:template>

</xsl:stylesheet>
