<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
		xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:loc="http://regis-web.systemsbiology.net/pepXML">
  <xsl:output method="text" />

  <xsl:template match="loc:mixture_model">
    <xsl:text>d </xsl:text>
    <xsl:value-of select ="@precursor_ion_charge"/><xsl:text> </xsl:text><xsl:value-of select ="@prior_probability"/><xsl:text>&#10;</xsl:text>
  </xsl:template>

  <xsl:template match="loc:spectrum_query/loc:search_result/loc:search_hit">
    <xsl:if test="@hit_rank='1'">
      <xsl:text>e </xsl:text><xsl:value-of select ="@peptide"/><xsl:text>&#10;</xsl:text>
      <xsl:text>c </xsl:text><xsl:value-of select ="../../@assumed_charge"/><xsl:text>&#10;</xsl:text>
      <xsl:text>r </xsl:text><xsl:value-of select ="@protein"/><xsl:text>&#10;</xsl:text>
      <xsl:for-each select="loc:alternative_protein">
	<xsl:text>r </xsl:text>
	<xsl:value-of select="@protein" /><xsl:text>&#10;</xsl:text>
      </xsl:for-each>
      <xsl:if test="loc:analysis_result/loc:peptideprophet_result/@probability!=''">
      <xsl:text>p </xsl:text><xsl:value-of select ="loc:analysis_result/loc:peptideprophet_result/@probability"/><xsl:text>&#10;</xsl:text>
      </xsl:if>
    </xsl:if>
  </xsl:template>

  <xsl:template match="text()"/>

</xsl:stylesheet>
